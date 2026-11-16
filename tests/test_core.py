import importlib.util
import struct
import tempfile
import hashlib
import xml.etree.ElementTree as ET
import unittest
from pathlib import Path
import zipfile
from unittest.mock import patch

MODULE=Path(__file__).resolve().parents[1]/'scripts'/'dfq.py'
spec=importlib.util.spec_from_file_location('dfq',MODULE)
dfq=importlib.util.module_from_spec(spec);spec.loader.exec_module(dfq)

class Core(unittest.TestCase):
    def test_generic_fragment_names_not_tied_to_vendor_order(self):
        self.assertEqual(dfq.fragment_id('G32159-1-10 PUC57'), '1-10')
        self.assertEqual(dfq.fragment_id('ProjectABC-2-27 sample'), '2-27')
        self.assertEqual(dfq.fragment_id('3-7'), '3-7')
        self.assertIsNone(dfq.fragment_id('unnumbered_fragment'))
    def test_parser_snapgene_reference_label_not_genotype(self):
        seq=b'A'*400+b'C'*400
        features=b'<Features><Feature name="1-C-T" type="misc_feature"><Segment range="500-500"/></Feature></Features>'
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'Sample-1-2 PUC57-KAN.dna'
            with p.open('wb') as f:
                for typ,val in [(0,b'\x00'+seq),(10,features)]:
                    f.write(bytes([typ])+struct.pack('>I',len(val))+val)
            rows=dfq.extract(p)
            self.assertEqual(len(rows),1)
            self.assertEqual(rows[0]['dna_base'],'C')
            dfq.apply_manifest(rows,{})
            self.assertEqual(rows[0]['reported_alt'],'T')
            self.assertEqual(rows[0]['dna_base'],'C')
    def test_multimap_gates(self):
        self.assertEqual(dfq.decide({'mapping':'MULTIMAP'})[0],'暂不判定')
        self.assertEqual(dfq.decide({'mapping':'MAPPED','evidence':'LOCAL_MAPPING_ONLY'})[0],'证据不足')
        self.assertEqual(dfq.decide({'mapping':'MAPPED','evidence':'UCSC_API_VERIFIED','mutation':'1-A-G','ccre_hits':[{'name':'abc'}]})[0],'优先评估修复')
    def test_non_synonymous_effect(self):
        # minus-strand, 1 exon; genomic C->T causes coding G->A: GTC -> ATC, V->I
        transcript={'name':'NM_test.1','name2':'GENE','strand':'-', 'cdsStart':100,'cdsEnd':106,
                    'exonStarts':'100,','exonEnds':'106,'}
        r={'chrom':'chrTest','hg38_pos':106,'ref_base':'C','mutation':'1-C-T','gene_hits':[transcript]}
        def fake_get(_endpoint,**_kw):
            # reverse complement = GTCATG, site zero 102 -> coding index 3: A in this case
            return {'dna':'CATGAC'}
        effect=dfq.coding_effect(r,fetch=fake_get)
        self.assertIsNotNone(effect)
        self.assertEqual(effect['effect'],'nonsynonymous')
    def test_cli_four_outputs_with_mocked_mapping_and_ucsc(self):
        seq=b'A'*400+b'C'*400
        xml=b'<Features><Feature name="1-C-T" type="misc_feature"><Segment range="500-500"/></Feature></Features>'
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'Sample-1-2 PUC57-KAN.dna'
            with p.open('wb') as out:
                for kind,contents in [(0,b'\x00'+seq),(10,xml)]:
                    out.write(bytes([kind])+struct.pack('>I',len(contents))+contents)
            def mapping(records,*args):
                records[0].update(mapping='MAPPED',chrom='chr17',hg38_pos=123456,
                                  strand='+',identity=100.0,top_hits=[{'chrom':'chr17','pos':123456}])
                return records
            def annotate(records,*args):
                records[0].update(evidence='UCSC_API_VERIFIED',ref_base='C',ref_context='ACA',
                                  repeat_hits=[],ccre_hits=[],gene_hits=[])
                return records
            with patch.object(dfq,'map_blast',mapping),patch.object(dfq,'get_annotations',annotate):
                self.assertEqual(dfq.main(['--input',d,'--db','unused','--out',str(Path(d)/'out'),'--online']),0)
            self.assertEqual({x.name for x in (Path(d)/'out').iterdir()},
                             {'review.xlsx','review.tsv','report.md','result.json'})
            self.assertIn('一般可暂不修复',(Path(d)/'out'/'report.md').read_text(encoding='utf-8'))

    def test_vendor_copy_last_step_preserves_original_and_never_overwrites(self):
        # Entirely synthetic vendor Excel, including a second worksheet to preserve.
        sheet = ('<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">'
                 '<dimension ref="A1:G2"/><sheetData><row r="1">'
                 '<c r="B1" t="inlineStr"><is><t>4K分段</t></is></c>'
                 '</row><row r="2">'
                 '<c r="B2" t="inlineStr"><is><t>Project-1-1</t></is></c>'
                 '<c r="E2" t="inlineStr"><is><t>1-A-G；2-C-T；3-G-A</t></is></c>'
                 '<c r="F2" t="inlineStr"><is><t>原结论</t></is></c>'
                 '</row></sheetData></worksheet>')
        other = b'<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main"><sheetData/></worksheet>'
        with tempfile.TemporaryDirectory() as d:
            source = Path(d) / 'vendor.xlsx'
            with zipfile.ZipFile(source,'w') as z:
                z.writestr('xl/worksheets/sheet1.xml',sheet)
                z.writestr('xl/worksheets/sheet2.xml',other)
            before = hashlib.sha256(source.read_bytes()).hexdigest()
            records = [
                dict(fragment='1-1',annotation='1-A-G',mutation='1-A-G',
                     decision='优先评估修复',reason='candidate',mapping='MAPPED',
                     evidence='UCSC_API_VERIFIED',chrom='chr17',hg38_pos=101,
                     ccre_hits=[{'name':'EH-test'}],gene_hits=[],repeat_hits=[]),
                dict(fragment='1-1',annotation='2-C-T',mutation='2-C-T',
                     decision='一般可暂不修复',reason='no annotation',mapping='MAPPED',
                     evidence='UCSC_API_VERIFIED',chrom='chr17',hg38_pos=102,
                     ccre_hits=[],gene_hits=[],repeat_hits=[],cpg={'delta':0,'reference_context':'ACA','alternate_context':'ATA'}),
            ]
            output = dfq.export_annotated_vendor_copy(source,Path(d),records)
            self.assertNotEqual(source,output)
            self.assertEqual(hashlib.sha256(source.read_bytes()).hexdigest(), before)
            with zipfile.ZipFile(output) as z:
                self.assertEqual(z.read('xl/worksheets/sheet2.xml'),other)
                ns={'m':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
                root=ET.fromstring(z.read('xl/worksheets/sheet1.xml'))
                cells={c.get('r'):''.join(c.itertext()) for c in root.findall('.//m:c',ns)}
                self.assertEqual(cells['E2'],'1-A-G；2-C-T；3-G-A')
                self.assertEqual(cells['F2'],'1修复；2不修复；3待核实')
                self.assertIn('候选调控元件EH-test',cells['G2'])
                self.assertIn('未找到对应的SnapGene',cells['G2'])
            first_hash = hashlib.sha256(output.read_bytes()).hexdigest()
            second=dfq.export_annotated_vendor_copy(source,Path(d),records)
            self.assertNotEqual(output,second)
            self.assertEqual(hashlib.sha256(output.read_bytes()).hexdigest(),first_hash)

    def test_final_calls_repair_nonsynonymous_and_withhold_uncertain(self):
        base=dict(fragment='1-2',annotation='2-C-T',mutation='2-C-T',mapping='MAPPED',
                  chrom='chrTest',hg38_pos=123456,evidence='UCSC_API_VERIFIED',
                  gene_hits=[],ccre_hits=[],repeat_hits=[],
                  decision='建议修复',reason='coding',cpg={'delta':0,'reference_context':'ACA','alternate_context':'ATA'})
        altered=dict(base,coding={'gene':'GENE_TEST','effect':'nonsynonymous',
                                  'transcript':'NM_TEST.1','aa':'A12V'})
        self.assertEqual(dfq.final_repair_call(altered)[0],'修复')
        self.assertIn('A12V',dfq.final_repair_call(altered)[1])
        unclear=dict(base,decision='核对原始注释',reason='BC1不等于报告1-A-T')
        self.assertEqual(dfq.final_repair_call(unclear)[0],'待核实')
        offline=dict(base,decision='证据不足',reason='UCSC未核验')
        self.assertEqual(dfq.final_repair_call(offline)[0],'待核实')
        repeat=dict(base,annotation='1-C-CACAC',mutation='1-C-CACAC',reported_ref='C',reported_alt='CACAC',
                    decision='条件修复',coding=None,repeat_hits=[{'repName':'(AC)n','repClass':'Simple_repeat'}])
        self.assertEqual(dfq.final_repair_call(repeat)[0],'不修复')
        self.assertIn('实验变量',dfq.final_repair_call(repeat)[1])
    def test_cpg_gain_loss_unaffected_and_minus_strand(self):
        samples=[
            ('1-A-G','CAA','+',1,'CGA','修复'),
            ('1-C-T','ACG','+',-1,'ATG','修复'),
            ('1-T-C','CAA','-',1,'CGA','修复'),
            ('1-T-C','CTC','+',0,'CCC','不修复'),
        ]
        for allele,triplet,strand,delta,expected,call in samples:
            with self.subTest(allele=allele,triplet=triplet,strand=strand):
                r=dict(fragment='1-1',annotation=allele,mutation=allele,mapping='MAPPED',
                       chrom='chrTest',hg38_pos=11,strand=strand,ref_base=triplet[1],
                       ref_context=triplet,evidence='UCSC_API_VERIFIED',
                       repeat_hits=[],ccre_hits=[],gene_hits=[])
                r['cpg']=dfq.cpg_change(r)
                self.assertIsNotNone(r['cpg'])
                self.assertEqual(r['cpg']['delta'],delta)
                self.assertEqual(r['cpg']['alternate_context'],expected)
                r['decision'],r['reason']=dfq.decide(r)
                self.assertEqual(dfq.final_repair_call(r)[0],call)
                if delta:self.assertIn('CpG',dfq.final_repair_call(r)[1])

    def test_cpg_reference_context_reuses_one_ucsc_sequence_request(self):
        r=dict(fragment='1-1',annotation='1-A-G',mutation='1-A-G',mapping='MAPPED',
               chrom='chrTest',hg38_pos=11,strand='+')
        calls=[]
        def fake_api(endpoint,**kwargs):
            calls.append((endpoint,kwargs))
            return {'dna':'CAA'} if endpoint=='getData/sequence' else {kwargs['track']:[]}
        with patch.object(dfq,'api',fake_api):
            dfq.get_annotations([r],online=True)
        self.assertEqual(r['ref_base'],'A')
        self.assertEqual(r['ref_context'],'CAA')
        self.assertEqual(r['evidence'],'UCSC_API_VERIFIED')
        self.assertEqual(sum(x[0]=='getData/sequence' for x in calls),1)
        self.assertEqual([(x[1]['start'],x[1]['end']) for x in calls if x[0]=='getData/sequence'],[(9,12)])
        self.assertEqual(dfq.cpg_change(r)['delta'],1)

    def test_cpg_missing_context_never_classified_as_no_change(self):
        r=dict(fragment='1-1',annotation='1-A-G',mutation='1-A-G',mapping='MAPPED',
               chrom='chrTest',hg38_pos=11,strand='+',ref_base='A',
               evidence='UCSC_API_VERIFIED',ccre_hits=[],repeat_hits=[],gene_hits=[])
        self.assertIsNone(dfq.cpg_change(r))
        r['cpg']=None
        r['decision'],r['reason']=dfq.decide(r)
        self.assertEqual(dfq.final_repair_call(r)[0],'待核实')
        r['ref_context']='CAA'
        r['annotation_warning']='BC1 is not mutation 1-A-G'
        self.assertIsNone(dfq.cpg_change(r))
        r['mutation']='1-A-AG'
        r.pop('annotation_warning')
        self.assertIsNone(dfq.cpg_change(r))

    def test_xlsx_opens_as_valid_package(self):
        with tempfile.TemporaryDirectory() as d:
            x=Path(d)/'result.xlsx'
            dfq.export_xlsx(x,[{'fragment':'F1','mutation':'1-A-G','decision':'建议修复'}])
            with zipfile.ZipFile(x) as z:
                self.assertIn('xl/worksheets/sheet1.xml',z.namelist())
                self.assertIn('建议修复',z.read('xl/worksheets/sheet1.xml').decode())

if __name__=='__main__':unittest.main()
