import importlib.util
import struct
import tempfile
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
                records[0].update(evidence='UCSC_API_VERIFIED',ref_base='C',
                                  repeat_hits=[],ccre_hits=[],gene_hits=[])
                return records
            with patch.object(dfq,'map_blast',mapping),patch.object(dfq,'get_annotations',annotate):
                self.assertEqual(dfq.main(['--input',d,'--db','unused','--out',str(Path(d)/'out'),'--online']),0)
            self.assertEqual({x.name for x in (Path(d)/'out').iterdir()},
                             {'review.xlsx','review.tsv','report.md','result.json'})
            self.assertIn('一般可暂不修复',(Path(d)/'out'/'report.md').read_text(encoding='utf-8'))
    def test_xlsx_opens_as_valid_package(self):
        with tempfile.TemporaryDirectory() as d:
            x=Path(d)/'result.xlsx'
            dfq.export_xlsx(x,[{'fragment':'F1','mutation':'1-A-G','decision':'建议修复'}])
            with zipfile.ZipFile(x) as z:
                self.assertIn('xl/worksheets/sheet1.xml',z.namelist())
                self.assertIn('建议修复',z.read('xl/worksheets/sheet1.xml').decode())

if __name__=='__main__':unittest.main()
