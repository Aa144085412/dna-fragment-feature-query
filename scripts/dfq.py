#!/usr/bin/env python3
"""DNA Fragment Feature Query: annotated variants -> hg38 BLAST -> UCSC -> evidence review.

Important: SnapGene sequence and reported mutant allele can be different; the
presence of a feature label is NOT evidence that sequencing confirmed an allele.
"""
from __future__ import annotations
import argparse
import csv
from collections import Counter, defaultdict
from datetime import datetime, timezone
import hashlib
import json
import re
import shutil
import struct
import subprocess
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
import zipfile
from pathlib import Path

DNA = set('ACGTN')
AL = re.compile(r'^([1-9]\d*)-([ACGT]+)-([ACGT]+)$')
ORD = re.compile(r'^([1-9]\d*)[-－]')
XNS = {'m':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
COLS = [('fragment','片段'),('mutation','报告突变'),('position','SnapGene坐标'),
 ('chrom','hg38染色体'),('hg38_pos','hg38位置'),('strand','方向'),
 ('mapping','BLAST定位状态'),('identity','比对一致性%'),('ref_base','hg38参考碱基'),
 ('reported_ref','报告REF'),('reported_alt','报告ALT'),('repeat','RepeatMasker'),
 ('ccre','ENCODE4 cCRE'),('genes','RefSeq基因'),('coding_effect','编码影响预测'),
 ('decision','修复建议'),('reason','判定依据'),('evidence','证据状态'),('ucsc','UCSC链接')]

def fragment_id(value):
    """Extract the final two numeric fragment segments, independent of order prefix.

    Examples: G32159-1-10 -> 1-10; Sample-2-8 -> 2-8; 3-7 -> 3-7.
    When there is no structured fragment ID, callers use the file stem.
    """
    candidates = re.findall(r'(?=(\d+-\d+))', str(value))
    return candidates[-1] if candidates else None


def snapgene(path):
    sequence,features = None,None
    with path.open('rb') as h:
        while (tag:=h.read(1)):
            size=h.read(4)
            if len(size)!=4: raise ValueError('Truncated SnapGene record')
            n=struct.unpack('>I',size)[0]
            if n>100_000_000: raise ValueError('SnapGene record too large')
            b=h.read(n)
            if len(b)!=n: raise ValueError('Incomplete SnapGene record')
            if tag==b'\x00': sequence=b[1:].decode('ascii').upper()
            if tag==b'\x0a': features=ET.fromstring(b)
    if not sequence or set(sequence)-DNA or features is None:
        raise ValueError(f'Missing DNA sequence or features: {path}')
    return sequence,features

def extract(path,flank=150,query_min=401):
    seq,root=snapgene(path)
    fragment=fragment_id(path.stem) or path.stem.split(' ')[0]
    sha=hashlib.sha256(path.read_bytes()).hexdigest()
    rows=[]
    for f in root.findall('.//Feature'):
        if f.get('type')!='misc_feature': continue
        label=f.get('name','')
        for q in f.findall('.//Q'):
            v=q.find('V')
            if q.get('name')=='label' and v is not None:label=v.get('text',label)
        if not ((ORD.match(label) and (label.count('-')>=2 or '缺失' in label)) or label=='BC1'):
            continue
        pos=f.find('.//Segment')
        if pos is None or not pos.get('range'):continue
        a,b=(int(v) for v in pos.get('range').split('-'))
        if not 1<=a<=b<=len(seq):raise ValueError(f'Invalid feature coordinates {a}-{b}')
        left=max(query_min,a-flank)
        if left>a:left=max(1,a-flank)
        right=min(len(seq),b+flank)
        rows.append(dict(fragment=fragment,annotation=label,mutation=label,position=a,end=b,
          sequence=seq[left-1:right],query_pos=a-left+1,dna_base=seq[a-1:b],
          source_file=path.name,source_sha256=sha))
    return rows

def manifest_from_xlsx(path):
    """Minimal Excel reader for vendor sheets: B=fragment, E=semicolon variant claims."""
    with zipfile.ZipFile(path) as z:
        strings=[]
        if 'xl/sharedStrings.xml' in z.namelist():
            strings=[''.join(n.itertext()) for n in ET.fromstring(z.read('xl/sharedStrings.xml')).findall('m:si',XNS)]
        sheets=sorted(k for k in z.namelist() if re.fullmatch(r'xl/worksheets/sheet\d+\.xml',k))
        if not sheets:raise ValueError('No .xlsx worksheet found')
        doc=ET.fromstring(z.read(sheets[0]))
    out={}
    for row in doc.findall('.//m:sheetData/m:row',XNS):
        cols={}
        for c in row.findall('m:c',XNS):
            col=re.match(r'[A-Z]+',c.get('r',''))
            v=c.find('m:v',XNS)
            inline=c.find('m:is',XNS)
            value=v.text if v is not None else ''.join(inline.itertext()) if inline is not None else ''
            value=value or ''
            if c.get('t')=='s':value=strings[int(value)]
            if col:cols[col[0]]=value
        name=fragment_id(cols.get('B',''))
        if name:
            claims={}
            for x in re.split('[;；]',cols.get('E','')):
                x=x.strip();m=ORD.match(x)
                if m:claims[int(m[1])]=x
            out[name]=claims
    return out

def apply_manifest(rows,manifest):
    for r in rows:
        idx=1 if r['annotation']=='BC1' else int(ORD.match(r['annotation'])[1])
        claim=manifest.get(r['fragment'],{}).get(idx)
        if claim:
            r['mutation']=claim
            r['manifest_claim']=claim
        if r['annotation']=='BC1':
            r['annotation_warning']='SnapGene only says BC1, while Excel describes variant 1; check raw trace'
        elif claim and r['annotation']!=claim and not ('缺失' in claim and '缺失' in r['annotation']):
            r['annotation_warning']=f"Excel: {claim}; SnapGene: {r['annotation']}"
        m=AL.match(r['mutation'])
        r['reported_ref']=m[2] if m else ''
        r['reported_alt']=m[3] if m else ''
    return rows

def map_blast(rows,blastn,db,threads=4):
    if not (shutil.which(blastn) or Path(blastn).is_file()):raise FileNotFoundError(f'blastn not found: {blastn}')
    fasta=''.join(f'>{i}\n{r["sequence"]}\n' for i,r in enumerate(rows))
    args=[blastn,'-db',db,'-query','-','-task','megablast','-dust','no',
          '-soft_masking','false','-perc_identity','90','-qcov_hsp_perc','70',
          '-max_target_seqs','25','-num_threads',str(threads),
          '-outfmt','6 qseqid sseqid pident length qstart qend sstart send bitscore qseq sseq']
    run=subprocess.run(args,input=fasta,text=True,capture_output=True,timeout=600)
    if run.returncode:raise RuntimeError('BLAST failed: '+run.stderr[:600])
    candidates=defaultdict(list)
    for line in run.stdout.splitlines():
        c=line.split('\t')
        if len(c)!=11 or not c[0].isdigit():continue
        i=int(c[0])
        if i>=len(rows):continue
        _,chrom,pid,length,qs,qe,ss,se,bits,qa,sa=c
        target=rows[i]['query_pos'];qindex=int(qs)-1;pos=int(ss)
        direction=1 if int(se)>=int(ss) else -1
        site=None
        for a,b in zip(qa,sa):
            if a!='-':qindex+=1
            if b!='-':
                if a!='-' and qindex==target:site=pos
                pos+=direction
        if site is None:continue
        coverage=(int(qe)-int(qs)+1)/len(rows[i]['sequence'])
        candidates[i].append(dict(chrom=chrom,pos=site,strand='+' if direction>0 else '-',
                         identity=float(pid),coverage=round(coverage,4),bits=float(bits),aligned=int(length)))
    for i,r in enumerate(rows):
        hits=sorted(candidates[i],key=lambda h:(h['bits'],h['coverage'],h['identity'],int('_' not in h['chrom'])),reverse=True)
        r['top_hits']=hits[:10]
        if not hits:r['mapping']='UNMAPPED';continue
        top=hits[0]
        tied=[h for h in hits[1:] if h['bits']>=top['bits']-1.0 and (h['chrom'],h['pos'])!=(top['chrom'],top['pos'])]
        r['mapping']='MULTIMAP' if tied or top['coverage']<.95 or top['identity']<98 else 'MAPPED'
        if tied and '_' not in top['chrom'] and all(h['chrom'].startswith(top['chrom']+'_') for h in tied):
            r['mapping']='MAPPED_WITH_ALT_CONTIG'
        r.update(chrom=top['chrom'],hg38_pos=top['pos'],strand=top['strand'],identity=top['identity'])
    return rows

def api(endpoint,**kwargs):
    url='https://api.genome.ucsc.edu/'+endpoint+'?'+urllib.parse.urlencode(kwargs)
    request=urllib.request.Request(url,headers={'User-Agent':'DNAFragmentFeatureQuery/0.1'})
    with urllib.request.urlopen(request,timeout=30) as r:return json.load(r)

def get_annotations(rows,online=False):
    for r in rows:
        if r.get('mapping') not in ('MAPPED','MAPPED_WITH_ALT_CONTIG'):
            r['evidence']='MAP_UNCERTAIN';continue
        if '_' in r['chrom']:
            r['evidence']='ALT_CONTIG; verify location manually';continue
        if not online:
            r['evidence']='LOCAL_MAPPING_ONLY; UCSC annotations not queried';continue
        chrom,pos=r['chrom'],r['hg38_pos']-1
        errors=[]
        try:r['ref_base']=api('getData/sequence',genome='hg38',chrom=chrom,start=pos,end=pos+1).get('dna','').upper()
        except Exception as e:errors.append('ref:'+str(e)[:120])
        for track,key in [('rmsk','repeat_hits'),('cCREregistry','ccre_hits'),('ncbiRefSeq','gene_hits')]:
            try:r[key]=api('getData/track',genome='hg38',track=track,chrom=chrom,start=pos,end=pos+1).get(track,[])
            except Exception as e:r[key]=None;errors.append(track+':'+str(e)[:120])
        if len(r.get('ref_base',''))!=1:errors.append('reference_base_not_verified')
        r['evidence']='UCSC_API_VERIFIED' if not errors else 'PARTIAL_UCSC: '+'; '.join(errors)
    return rows

CODE={c:a for c,a in zip(
 ['TTT','TTC','TTA','TTG','TCT','TCC','TCA','TCG','TAT','TAC','TAA','TAG','TGT','TGC','TGA','TGG',
  'CTT','CTC','CTA','CTG','CCT','CCC','CCA','CCG','CAT','CAC','CAA','CAG','CGT','CGC','CGA','CGG',
  'ATT','ATC','ATA','ATG','ACT','ACC','ACA','ACG','AAT','AAC','AAA','AAG','AGT','AGC','AGA','AGG',
  'GTT','GTC','GTA','GTG','GCT','GCC','GCA','GCG','GAT','GAC','GAA','GAG','GGT','GGC','GGA','GGG'],
 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG')}
COMP=str.maketrans('ACGT','TGCA')

def coding_effect(r,fetch=api):
    m=AL.match(r['mutation'])
    if not m or len(m[2])!=1 or len(m[3])!=1:return None
    genomic_ref=m[2].translate(COMP) if r.get('strand')=='-' else m[2]
    genomic_alt=m[3].translate(COMP) if r.get('strand')=='-' else m[3]
    if genomic_ref!=r.get('ref_base'):return None
    p=r['hg38_pos']-1
    for t in r.get('gene_hits') or []:
        if not t.get('name','').startswith('NM_'):continue
        a=t.get('exonStarts',[]);b=t.get('exonEnds',[])
        if isinstance(a,str):a=[int(x) for x in a.strip(',').split(',') if x]
        if isinstance(b,str):b=[int(x) for x in b.strip(',').split(',') if x]
        cds0,cds1=int(t.get('cdsStart',0)),int(t.get('cdsEnd',0))
        intervals=[(max(int(x),cds0),min(int(y),cds1)) for x,y in zip(a,b)]
        intervals=[(x,y) for x,y in intervals if x<y]
        reverse=t.get('strand')=='-'
        if reverse:intervals.reverse()
        if not any(x<=p<y for x,y in intervals):continue
        cds=[];offset=None;length=0
        for x,y in intervals:
            seq=fetch('getData/sequence',genome='hg38',chrom=r['chrom'],start=x,end=y).get('dna','').upper()
            if len(seq)!=y-x:raise ValueError('UCSC returned truncated CDS')
            if reverse:seq=seq.translate(COMP)[::-1]
            if x<=p<y:offset=length+(y-1-p if reverse else p-x)
            cds.append(seq);length+=len(seq)
        seq=''.join(cds)
        if offset is None:continue
        ref_codon=seq[offset//3*3:offset//3*3+3]
        if len(ref_codon)!=3:return None
        alt=genomic_alt.translate(COMP) if reverse else genomic_alt
        alt_codon=ref_codon[:offset%3]+alt+ref_codon[offset%3+1:]
        if ref_codon not in CODE or alt_codon not in CODE:return None
        aa=CODE[ref_codon]+str(offset//3+1)+CODE[alt_codon]
        return dict(gene=t.get('name2',''),transcript=t['name'],codon=ref_codon+'>'+alt_codon,
                    aa=aa,effect='synonymous' if CODE[ref_codon]==CODE[alt_codon] else 'nonsynonymous')
    return None

def decide(r):
    if r.get('mapping') not in ('MAPPED','MAPPED_WITH_ALT_CONTIG'):
        return '暂不判定','未可靠定位或全基因组多重命中'
    if r.get('annotation_warning'):
        return '核对原始注释',r['annotation_warning']
    genomic_ref=r.get('reported_ref','').translate(COMP) if r.get('strand')=='-' else r.get('reported_ref','')
    if len(genomic_ref)==1 and r.get('ref_base') and genomic_ref!=r['ref_base']:
        return '核对原始注释',f"报告参考碱基（换算至hg38正链）{genomic_ref} 与hg38 {r['ref_base']} 不符"
    if r.get('evidence')!='UCSC_API_VERIFIED':
        return '证据不足','UCSC功能注释未完整查询，不能把缺失的注释当作阴性'
    if r.get('coding_error'):
        return '证据不足','转录本CDS查询或翻译出错: '+str(r['coding_error'])
    effect=r.get('coding')
    if effect and effect['effect']=='nonsynonymous':
        return '建议修复',f"编码非同义变化预测 {effect['transcript']}: {effect['aa']}"
    if r.get('ccre_hits'):
        return '优先评估修复','与ENCODE4候选cCRE重叠；功能未被此查询证明'
    if '缺失' in r['mutation'] or (r.get('reported_ref') and len(r['reported_ref'])!=len(r.get('reported_alt',''))):
        return '条件修复','插入/缺失改变片段长度，重复区域需归一化并核对测序峰图'
    if effect and effect['effect']=='synonymous':
        return '一般可暂不修复','同义编码替换预测；仍须考虑实验设计'
    for g in r.get('gene_hits') or []:
        a=g.get('exonStarts',[]);b=g.get('exonEnds',[])
        if isinstance(a,str):a=[int(x) for x in a.strip(',').split(',') if x]
        if isinstance(b,str):b=[int(x) for x in b.strip(',').split(',') if x]
        if any(x<=r['hg38_pos']-1<y for x,y in zip(a,b)):
            return '需人工评估','位于RefSeq外显子，但尚无可解释的编码变化'
    return '一般可暂不修复','未查询到该位点的编码外显子/cCRE重叠；不能推出无功能'

def export_xlsx(path,records):
    """Minimal valid XLSX, standard-library only, no private source data included."""
    from xml.sax.saxutils import escape
    def column(n):
        s=''
        while n:n,mod=divmod(n-1,26);s=chr(65+mod)+s
        return s
    rows=[[v for _,v in COLS]]+[[str(r.get(k,'')) for k,_ in COLS] for r in records]
    data=[]
    for i,row in enumerate(rows,1):
        cells=[]
        for j,v in enumerate(row,1):
            if v.startswith(('=','+','-','@')):v="'"+v
            cells.append(f'<c r="{column(j)}{i}" t="inlineStr"><is><t>{escape(v)}</t></is></c>')
        data.append(f'<row r="{i}">'+''.join(cells)+'</row>')
    xml='<?xml version="1.0" encoding="UTF-8"?><worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main"><sheetViews><sheetView workbookViewId="0"><pane ySplit="1" state="frozen" topLeftCell="A2" activePane="bottomLeft"/></sheetView></sheetViews><sheetData>'+''.join(data)+'</sheetData><autoFilter ref="A1:S'+str(len(rows))+'"/></worksheet>'
    with zipfile.ZipFile(path,'w',zipfile.ZIP_DEFLATED) as z:
        z.writestr('[Content_Types].xml','<?xml version="1.0"?><Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types"><Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/><Default Extension="xml" ContentType="application/xml"/><Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/><Override PartName="/xl/worksheets/sheet1.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/></Types>')
        z.writestr('_rels/.rels','<?xml version="1.0"?><Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="xl/workbook.xml"/></Relationships>')
        z.writestr('xl/workbook.xml','<?xml version="1.0"?><workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"><sheets><sheet name="hg38逐条审查" sheetId="1" r:id="rId1"/></sheets></workbook>')
        z.writestr('xl/_rels/workbook.xml.rels','<?xml version="1.0"?><Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" Target="worksheets/sheet1.xml"/></Relationships>')
        z.writestr('xl/worksheets/sheet1.xml',xml)


def final_repair_call(r):
    """Binary vendor recommendation after the existing scientific review gates.

    Unresolved evidence is never silently converted into 'do not repair'.
    """
    prior, note = r['decision'], r['reason']
    position = f"{r['chrom']}:{r['hg38_pos']:,}" if r.get('chrom') and r.get('hg38_pos') else 'hg38坐标未确认'
    prefix = position + '，'
    if prior in ('暂不判定', '核对原始注释', '证据不足'):
        return '待核实', prefix + note
    effect = r.get('coding') or {}
    if effect.get('effect') == 'nonsynonymous':
        return '修复', prefix + f"{effect.get('gene','')}编码区，按{effect['transcript']}预测为{effect['aa']}非同义替换；建议恢复设计序列（预测不等于功能验证）"
    ccre = r.get('ccre_hits') or []
    if ccre:
        ids = ', '.join(x.get('name','未知') for x in ccre)
        return '修复', prefix + f'重叠ENCODE4候选调控元件{ids}；为保留潜在调控序列建议修复，但尚未证实该突变改变功能'
    repeats = r.get('repeat_hits') or []
    rep = ', '.join(dict.fromkeys(x.get('repName','?') for x in repeats))
    is_indel = ('缺失' in r['mutation'] or
                bool(r.get('reported_ref')) and len(r['reported_ref']) != len(r.get('reported_alt','')))
    if is_indel:
        if repeats and all(x.get('repClass') in ('SINE', 'Simple_repeat') for x in repeats):
            return '不修复', prefix + f'位于{rep}重复序列，改变局部重复长度；未查询到编码损害或cCRE重叠。仅在重复长度不是实验变量时可不修复，原始测序及indel归一化仍需核对'
        return '待核实', prefix + '插入/缺失会改变序列长度，当前证据无法排除结构或功能影响；先核对测序及实验设计'
    if effect.get('effect') == 'synonymous':
        return '不修复', prefix + f"{effect.get('gene','')}编码区预测为同义替换；未见蛋白氨基酸改变，但不能排除RNA或调控影响"
    pos0 = r['hg38_pos'] - 1
    in_exon, coding_exon, genes = False, False, set()
    for g in r.get('gene_hits') or []:
        genes.add(g.get('name2') or g.get('name') or '未知基因')
        starts, ends = g.get('exonStarts', []), g.get('exonEnds', [])
        if isinstance(starts, str):starts = [int(x) for x in starts.strip(',').split(',') if x]
        if isinstance(ends, str):ends = [int(x) for x in ends.strip(',').split(',') if x]
        if any(a <= pos0 < b for a,b in zip(starts, ends)):
            in_exon = True
            if int(g.get('cdsStart', 0)) <= pos0 < int(g.get('cdsEnd', 0)):
                coding_exon = True
    if coding_exon:
        return '待核实', prefix + '位于蛋白编码外显子，但缺少可靠的密码子影响解释；不得直接判为无功能'
    if in_exon:
        return '不修复', prefix + f"位于{', '.join(sorted(genes))}非编码外显子；未见蛋白编码改变或已查询cCRE重叠。RNA功能影响尚未排除，若研究RNA通路需重新评估"
    if repeats:
        location = f'位于{rep}重复元件内'
    elif genes:
        location = f"位于{', '.join(sorted(genes))}转录区域的非外显子部分"
    else:
        location = '未发现已查询RefSeq转录本或RepeatMasker重叠'
    return '不修复', prefix + location + '；未见已查询编码外显子或ENCODE4 cCRE重叠，暂无明确损害证据（不能证明没有功能）'


def export_annotated_vendor_copy(source, output_dir, rows):
    """Last step: produce a fresh XLSX with original sheet/layout plus F/G review.

    Leave the input untouched. Preserve all other zip members and worksheets;
    refuse to overwrite even an earlier generated result.
    """
    source, output_dir = Path(source), Path(output_dir)
    outcome = output_dir / (source.stem + '_hg38修复说明.xlsx')
    number = 2
    while outcome.exists() or outcome.resolve() == source.resolve():
        outcome = output_dir / (source.stem + f'_hg38修复说明_{number}.xlsx')
        number += 1
    lookup = {}
    for r in rows:
        m = ORD.match(r['annotation'])
        index = 1 if r['annotation'] == 'BC1' else int(m[1])
        key = (r['fragment'], index)
        if key in lookup:raise ValueError(f'Duplicate annotated mutation: {key}')
        lookup[key] = r
    with zipfile.ZipFile(source) as original:
        names = original.namelist()
        worksheets = sorted(x for x in names if re.fullmatch(r'xl/worksheets/sheet\d+\.xml',x))
        if not worksheets:raise ValueError('XLSX contains no worksheet')
        target = worksheets[0]
        strings = []
        if 'xl/sharedStrings.xml' in names:
            strings = [''.join(n.itertext()) for n in ET.fromstring(original.read('xl/sharedStrings.xml')).findall('m:si',XNS)]
        doc = ET.fromstring(original.read(target))
        body = doc.find('m:sheetData',XNS)
        if body is None:raise ValueError('XLSX has no sheetData')
        def read_cell(row, column):
            cell = next((c for c in row.findall('m:c',XNS) if re.match('[A-Z]+',c.get('r','')).group() == column),None)
            if cell is None:return ''
            value = cell.find('m:v',XNS)
            inline = cell.find('m:is',XNS)
            raw = value.text if value is not None else ''.join(inline.itertext()) if inline is not None else ''
            return strings[int(raw)] if cell.get('t') == 's' and raw else (raw or '')
        def put(row, column, value):
            address = column + row.get('r')
            cell = next((c for c in row.findall('m:c',XNS) if c.get('r') == address), None)
            if cell is None:
                cell = ET.Element('{' + XNS['m'] + '}c', {'r':address})
                siblings = list(row)
                order = lambda x: sum((ord(c)-64)*26**i for i,c in enumerate(re.match('[A-Z]+',x.get('r','A')).group()[::-1]))
                spot = next((i for i,x in enumerate(siblings) if x.tag.endswith('}c') and order(x)>order(cell)),len(siblings))
                row.insert(spot,cell)
            for child in list(cell):cell.remove(child)
            cell.set('t','inlineStr')
            if value.startswith(('=','+','-','@')):value="'"+value
            ET.SubElement(ET.SubElement(cell,'{'+XNS['m']+'}is'),'{'+XNS['m']+'}t').text=value
        updated = 0
        for excel_row in body.findall('m:row', XNS):
            number = int(excel_row.get('r','0'))
            if number == 1:
                put(excel_row,'F','修复结论')
                put(excel_row,'G','每处突变的修复／不修复原因（hg38）')
                continue
            fragment = fragment_id(read_cell(excel_row,'B'))
            if not fragment:continue
            claims = [v.strip() for v in re.split('[;；]',read_cell(excel_row,'E')) if v.strip()]
            decisions, explanations = [], []
            for claim in claims:
                m = ORD.match(claim)
                if not m:continue
                index = int(m[1]); key=(fragment,index)
                record = lookup.get(key)
                if record is None:
                    call,reason='待核实','未找到对应的SnapGene突变注释，不能判断是否影响功能'
                else:
                    call,reason=final_repair_call(record)
                decisions.append(f'{index}{call}')
                explanations.append(f'{claim}｜{call}：{reason}')
            if not decisions:continue
            put(excel_row,'F','；'.join(decisions))
            put(excel_row,'G','\n'.join(explanations))
            updated += 1
        if not updated:raise ValueError('No vendor Excel rows could be matched; review not exported')
        dimension=doc.find('m:dimension',XNS)
        if dimension is not None:
            ref=dimension.get('ref','')
            if ':' in ref:
                beginning,end=ref.split(':',1)
                m=re.match(r'[A-Z]+(\d+)',end)
                if m and end[0] < 'G':dimension.set('ref',beginning+':G'+m[1])
        ET.register_namespace('',XNS['m'])
        ET.register_namespace('r','http://schemas.openxmlformats.org/officeDocument/2006/relationships')
        edited=ET.tostring(doc,encoding='utf-8',xml_declaration=True)
        with zipfile.ZipFile(outcome,'x') as dest:
            for info in original.infolist():
                dest.writestr(info,edited if info.filename == target else original.read(info.filename))
    return outcome

def main(argv=None):
    ap=argparse.ArgumentParser(description='DNA Fragment Feature Query hg38 review')
    ap.add_argument('--input',required=True,help='Directory of mutation-annotated SnapGene .dna files')
    ap.add_argument('--xlsx',help='Optional vendor XLSX, fragment in B and mutation list in E')
    ap.add_argument('--blastn',default='blastn',help='Existing NCBI BLAST+ executable')
    ap.add_argument('--db',required=True,help='Existing indexed hg38 BLAST database prefix')
    ap.add_argument('--out',default='dfq-results')
    ap.add_argument('--threads',type=int,default=4)
    ap.add_argument('--flank',type=int,default=150)
    ap.add_argument('--query-min',type=int,default=1,help='First genomic-insert base in plasmid (default 1; use 401 for pUC inserts starting at 401)')
    ap.add_argument('--online',action='store_true',help='Allow UCSC REST requests for annotations')
    args=ap.parse_args(argv)
    paths=sorted(Path(args.input).glob('*.dna'))
    if not paths:ap.error('--input contains no .dna files')
    rows=[]
    for f in paths:rows+=extract(f,args.flank,args.query_min)
    if not rows:ap.error('No mutation-labelled features found')
    claims=manifest_from_xlsx(Path(args.xlsx)) if args.xlsx else {}
    apply_manifest(rows,claims)
    observed={(r['fragment'],1 if r['annotation']=='BC1' else int(ORD.match(r['annotation'])[1])) for r in rows}
    unmatched=[{'fragment':fragment,'mutation':mutation} for fragment,entries in claims.items()
              for index,mutation in entries.items() if (fragment,index) not in observed]
    map_blast(rows,args.blastn,args.db,args.threads)
    get_annotations(rows,args.online)
    for r in rows:
        if r.get('mapping') in ('MAPPED','MAPPED_WITH_ALT_CONTIG'):
            chrom,pos=r['chrom'],r['hg38_pos']
            r['ucsc']=f'https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position={chrom}%3A{max(1,pos-100)}-{pos+100}'
        if r.get('evidence')=='UCSC_API_VERIFIED' and r.get('gene_hits'):
            try:r['coding']=coding_effect(r)
            except Exception as exc:r['coding_error']=str(exc)
        r['decision'],r['reason']=decide(r)
        r['repeat']=', '.join(x.get('repName','?')+' ('+x.get('repClass','?')+')' for x in r.get('repeat_hits') or [])
        r['ccre']=', '.join(x.get('name','?')+' ('+x.get('cCRE_class','?')+')' for x in r.get('ccre_hits') or [])
        r['genes']=', '.join(sorted({x.get('name2',x.get('name','?')) for x in r.get('gene_hits') or []}))
        r['coding_effect']=(r.get('coding') or {}).get('aa','')
    out=Path(args.out);out.mkdir(parents=True,exist_ok=True)
    report=dict(schema='dfq/0.1',reference='hg38',timestamp_utc=datetime.now(timezone.utc).isoformat(),
         truth='reported annotation only; no original sequencing trace validation',count=len(rows),
         source_files=len(paths),online=args.online,blast_database=args.db,unmatched_manifest_claims=unmatched,records=rows)
    (out/'result.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
    with (out/'review.tsv').open('w',encoding='utf-8-sig',newline='') as h:
        w=csv.writer(h,delimiter='\t');w.writerow([c for _,c in COLS]);w.writerows([[r.get(k,'') for k,_ in COLS] for r in rows])
    export_xlsx(out/'review.xlsx',rows)
    counter=Counter(r['decision'] for r in rows)
    md=['# DNA Fragment Feature Query — hg38逐条审查','',
        f'- 片段数 {len(paths)}；突变记录 {len(rows)}；基因组组装 hg38',
        '- **注意：突变来自Excel/注释，不是由本工具独立验证的实测等位基因。**',
        '- 不把cCRE当成已证实功能，也不把无注释当成无功能。',
        f'- Excel 中有 {len(unmatched)} 条未匹配的突变报告，请在 result.json 查看。','',
        '## 结论汇总','']+[f'- {k}: {v}' for k,v in counter.items()]+[
        '', '## 逐条结果','', '|片段|标签|坐标|重复|cCRE|编码影响|建议|','|---|---|---|---|---|---|---|']
    for r in rows:
        values=[r['fragment'],r['mutation'],f"{r.get('chrom','?')}:{r.get('hg38_pos','?')}",
                r['repeat'],r['ccre'],r['coding_effect'],r['decision']]
        md.append('|'+ '|'.join(str(v).replace('|','/')[:120] for v in values)+'|')
    md+=['','## 数据来源','', '本地 NCBI BLAST+ hg38 索引，UCSC /getData/sequence /getData/track。详见 result.json。']
    (out/'report.md').write_text('\n'.join(md)+'\n',encoding='utf-8')
    vendor_copy = export_annotated_vendor_copy(Path(args.xlsx),out,rows) if args.xlsx else None
    print(json.dumps({'files':len(paths),'variants':len(rows),'decisions':counter,'output':str(out.resolve()),
                       'vendor_review_copy':str(vendor_copy) if vendor_copy else None,
                       'online':args.online},ensure_ascii=False))
    return 0

if __name__=='__main__':raise SystemExit(main())
