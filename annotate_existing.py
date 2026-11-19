#!/usr/bin/env python3
"""Annotate retained Tortuga candidates; optionally repeat BOLD validation, never assembly."""
import argparse
from datetime import datetime, timezone
import fcntl
import json
import re
from pathlib import Path
import shutil
from types import SimpleNamespace

import mito_pipeline as m
import compact_outputs as c


def saved_options(config, cli):
    args=m.parser().parse_args([])
    for key,value in config.get('options',{}).items():
        if hasattr(args,key):
            setattr(args,key,value)
    for key in ['mitos_refdir','mitos_refseqver','bold_fasta','bold_metadata','threads','annotation_timeout','bold_max_hits','bold_search_max_hits']:
        value=getattr(cli,key,None)
        if value is not None:
            setattr(args,key,value)
    for key in ['mitos_refdir','bold_fasta','bold_metadata','env_root','conda','itsme_script']:
        if getattr(args,key,None):
            setattr(args,key,Path(getattr(args,key)))
    args.resume=True
    args.dry_run=False
    args.annotation_policy='best'
    # Existing commands are the actual per-tool environments used by the run.
    commands=config.get('commands',{})
    if cli.config:
        commands=m.resolve_commands(SimpleNamespace(**{**vars(args),'config':cli.config}),['mitos2','blastn','makeblastdb'])[0]
    elif not commands.get('mitos2'):
        commands=m.resolve_commands(args,['mitos2','blastn','makeblastdb'])[0]
    for tool in ['mitos2','blastn','makeblastdb']:
        if commands.get(tool) and (not isinstance(commands[tool],list) or not all(isinstance(x,str) for x in commands[tool])):
            raise ValueError(f'Invalid saved command for {tool}')
    if args.bold_max_hits<1 or args.bold_search_max_hits<args.bold_max_hits:
        raise ValueError('Invalid saved BOLD search limits')
    return args,commands


def retained_sequence(sample):
    for path in [sample/'MitoLooper/best_sequence.fasta',sample/'final/best_sequence.fasta']:
        if path.is_file() and path.stat().st_size:
            try:
                path.resolve().relative_to(sample.resolve())
            except ValueError:
                raise ValueError('Retained candidate points outside its sample directory')
            return path
    return None


def annotation_signature(source,code,args,commands):
    return {'input_sha256':m.digest_file(source),'mt_code':code,
            'mitos_refdir':str(args.mitos_refdir),'mitos_refseqver':args.mitos_refseqver,
            'command':commands.get('mitos2',[])}


def publish_annotation(sample,result,keep):
    """Export owned MITOS outputs before removing only this attempt's scratch."""
    if not result.get('directory'):
        return
    output=Path(result['directory'])
    mapping={}
    for path in c.owned_files(output):
        category=c.kind(path)
        if path.is_file() and category and not path.is_symlink():
            target=c.export_file(path,sample/'logs' if category=='log' else sample/'MitoLooper','mitos')
            mapping[str(path)]=str(target)
            if target.suffix.lower() in {'.gff','.gff3'}:
                c.annotation_table(target)
    annotation=sample/'MitoLooper/annotation_summary.json'
    if annotation.is_file():
        c.atomic_json(annotation,c.rewrite(c.read_json(annotation),mapping))
    if result['ok'] and not keep:
        # This is a freshly owned attempt directory, never a raw-read/reference path.
        work=output.parent
        work.resolve().relative_to((sample/'05_mitos2').resolve())
        shutil.rmtree(work)


def archive_acceptance(sample,history):
    mapping={}
    for name in ['accepted_candidate.fasta','circular_target_candidate.fasta']:
        for path in [sample/'MitoLooper'/name,sample/'final'/name]:
            if path.is_file() and not path.is_symlink():
                target=c.export_file(path,history,'previous')
                mapping[str(path)]=str(target)
                path.unlink()
    manifest=sample/'final/outputs_manifest.tsv'
    if manifest.is_file() and mapping:
        import csv
        with manifest.open(newline='') as handle:
            reader=csv.DictReader(handle,delimiter='\t')
            fields=list(reader.fieldnames or [])
            records=list(reader)
        for row in records:
            row['retained_file']=mapping.get(row['retained_file'],row['retained_file'])
        c.table(manifest,records,fields)


def recovered_status(summary):
    if summary.get('accepted_by_policy'):
        strict=summary.get('strict_complete_candidate')
        if summary.get('annotation_status')=='skipped_unknown_code':
            return 'accepted_by_policy_annotation_unavailable'
        return ('recovered_annotated' if strict else 'accepted_by_policy_annotated') if summary.get('annotation_ok') else ('recovered_annotation_failed' if strict else 'accepted_by_policy_annotation_failed')
    return 'partial_target_confirmed' if summary.get('target_status')=='confirmed_target' else 'unresolved_candidate'


def postprocess(root,cli):
    root=Path(root)
    from Bio import SeqIO
    from Bio.Data import CodonTable
    rows=[]
    metadata=None
    db=None
    db_stamp=None
    for path in sorted(root.glob('*/final/summary.json')):
        sample=path.parent.parent
        if cli.only_sample and sample.name!=cli.only_sample:
            continue
        if sample.is_symlink() or (sample/'MitoLooper').is_symlink() or not (sample/'run_config.json').is_file() or not c.finished_sample(sample):
            continue
        summary=c.read_json(path)
        if summary.get('pipeline')=='itsme':
            continue
        source=retained_sequence(sample)
        if source is None:
            rows.append({'sample':sample.name,'status':'skipped_no_candidate'})
            continue
        print(f'\n=== {sample.name}: saved mitochondrial candidate ===',flush=True)
        try:
            for directory in ['final','MitoLooper','05_mitos2','postprocess_history','candidates','logs','state']:
                if (sample/directory).is_symlink():
                    raise ValueError(f'Refusing linked output directory: {directory}')
            config=c.read_json(sample/'run_config.json')
            args,commands=saved_options(config,cli)
            candidate=dict(summary.get('best_candidate') or {})
            candidate.update(fasta=str(source),candidate_id=candidate.get('candidate_id') or 'retained_candidate',
                             circular=bool(summary.get('circular')),accepted=bool(summary.get('accepted_by_policy')))
            if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*',candidate['candidate_id']):
                raise ValueError('Unsafe saved candidate ID')
            records=list(SeqIO.parse(source,'fasta'))
            if len(records)!=1 or not records[0].seq:
                raise ValueError('Annotation requires one nonempty retained candidate sequence')
            saved_sample=config.get('sample')
            saved_code=saved_sample.get('mt_code') if isinstance(saved_sample,dict) else None
            code=cli.mt_code if cli.mt_code is not None else summary.get('mt_code') or saved_code
            code=int(code) if code is not None else None
            if code is not None and int(code) not in CodonTable.unambiguous_dna_by_id:
                raise ValueError('Invalid mitochondrial genetic code')
            if args.threads<1 or args.annotation_timeout<=0:
                raise ValueError('Threads and annotation timeout must be positive')
            if not args.mitos_refdir or not args.mitos_refseqver or not (args.mitos_refdir/args.mitos_refseqver).is_dir():
                raise ValueError('MITOS reference release unavailable; supply --mitos-refdir and --mitos-refseqver')
            runner=m.Runner(sample,commands,resume=False)
            signature=annotation_signature(source,code,args,commands)
            cache=sample/'MitoLooper/postprocess_annotation.json'
            old=c.read_json(cache) if cache.is_file() else {}
            annotated=sample/'MitoLooper/annotation_summary.json'
            published_gff=c.read_json(annotated).get('gff','') if annotated.is_file() else ''
            cached=(not cli.force and old.get('signature')==signature and old.get('status')=='completed'
                    and bool(published_gff) and (Path(published_gff).is_file() or (sample/'MitoLooper'/Path(published_gff).name).is_file()))
            if cached and not cli.recheck_bold:
                print('  MITOS: cached published annotation',flush=True)
                rows.append({'sample':sample.name,'status':'cached'})
                continue
            history=sample/'postprocess_history'
            history.mkdir(exist_ok=True)
            attempt=history/f'attempt_{len(list(history.glob("attempt_*")))+1:03d}'
            attempt.mkdir()
            c.atomic_json(attempt/'summary_before.json',summary)
            if (sample/'stages.json').is_file():
                old_stages=c.read_json(sample/'stages.json')
            else:
                import csv
                with (sample/'stages.tsv').open(newline='') as handle:
                    old_stages=list(csv.DictReader(handle,delimiter='\t'))
            c.atomic_json(attempt/'stages_before.json',{'stages':old_stages})
            recheck_error=None
            if cli.recheck_bold:
                try:
                    if not args.bold_fasta or not args.bold_metadata or not args.bold_fasta.is_file() or not args.bold_metadata.is_file():
                        raise ValueError('BOLD input files unavailable; supply --bold-fasta and --bold-metadata')
                    stamp={'fasta':m.stamp(args.bold_fasta),'metadata':m.stamp(args.bold_metadata)}
                    if db_stamp is not None and db_stamp!=stamp:
                        raise ValueError('Selected samples use different BOLD datasets; process separately')
                    if metadata is None:
                        loaded_metadata=m.load_bold_metadata(args.bold_metadata)
                        dbdir=root/'postprocess_bold_reference'
                        if dbdir.is_symlink():
                            raise ValueError('Refusing linked postprocess BOLD directory')
                        dbdir.mkdir(exist_ok=True)
                        db=dbdir/'coi'
                        brunner=m.Runner(root,commands,resume=True)
                        if not brunner.run('postprocess_bold_index','makeblastdb',['-in',args.bold_fasta,'-dbtype','nucl','-parse_seqids','-out',db],dbdir,args.validation_timeout*60,[dbdir/'coi.nsq']):
                            raise ValueError('BOLD index preparation failed')
                        metadata=loaded_metadata
                        db_stamp=stamp
                        c.atomic_json(dbdir/'source.json',{**stamp,'metadata_sha256':m.digest_file(args.bold_metadata)})
                    query=source
                    query_source='retained_candidate_against_COI_database'
                    protein_hits=sample/'candidates'/candidate['candidate_id']/'proteins.tsv'
                    if protein_hits.is_file():
                        hits=[row for row in m.blast_rows(protein_hits) if m.normalize_gene(row['qseqid'])=='cox1']
                        if hits:
                            hit=max(hits,key=lambda row:float(row['bitscore']))
                            start,end=int(hit['sstart']),int(hit['send'])
                            sequence=records[0].seq+records[0].seq if candidate['circular'] else records[0].seq
                            if 1<=min(start,end)<=max(start,end)<=len(sequence):
                                coi=sequence[min(start,end)-1:max(start,end)]
                                if start>end:
                                    coi=coi.reverse_complement()
                                if len(coi)>=args.coi_min_bp:
                                    query=attempt/'coi_query.fasta'
                                    m.fasta_write(query,'COI',coi)
                                    query_source='saved_protein_localization'
                    decision=m.screen_bold(runner,'postprocess_bold',query,db,attempt,metadata,summary.get('target_taxon',''),args)
                    decision['query_source']=query_source
                    c.atomic_json(attempt/'bold_evidence.json',decision)
                    if decision.get('validation_error'):
                        raise ValueError(decision['validation_error'])
                    summary['bold_recheck_status']='completed'
                    summary['bold_recheck_reason']=decision['reason']
                    candidate['target_status']=decision['status']
                    candidate['target_reason']=decision['reason']
                    if candidate.get('validation_error')=='COI/BOLD screening failed':
                        candidate.pop('validation_error',None)
                    m.apply_acceptance(candidate,args)
                    summary.update(target_status=decision['status'],target_reason=decision['reason'],
                                   accepted_by_policy=candidate['accepted'],strict_complete_candidate=candidate['strict_complete_candidate'],
                                   winning_method=candidate.get('method','') if candidate['accepted'] else '')
                    evidence=sample/'candidates'/candidate['candidate_id']
                    evidence.mkdir(parents=True,exist_ok=True)
                    for name in ['bold_evidence.json','bold_hits.tsv']:
                        if (evidence/name).is_file():
                            shutil.copyfile(evidence/name,attempt/('previous_'+name))
                        shutil.copyfile(attempt/name,evidence/name)
                    c.atomic_json(evidence/'assessment.json',candidate)
                    if not candidate['accepted']:
                        archive_acceptance(sample,attempt)
                    else:
                        shutil.copyfile(source,sample/'MitoLooper/accepted_candidate.fasta')
                        if candidate['strict_complete_candidate']:
                            shutil.copyfile(source,sample/'MitoLooper/circular_target_candidate.fasta')
                except (ValueError,OSError,RuntimeError) as error:
                    recheck_error=str(error)
                    summary['bold_recheck_status']='failed'
                    summary['bold_recheck_reason']=recheck_error
                    print(f'  BOLD recheck failed: {error}; continuing MITOS annotation',flush=True)
            if cached:
                result={'ok':True,'status':'completed','scope':'accepted_candidate' if candidate['accepted'] else 'exploratory_candidate',
                        'candidate_id':candidate['candidate_id'],'directory':''}
            else:
                result=m.annotate_mito_candidate(candidate,code,sample,runner,args)
                publish_annotation(sample,result,cli.keep_intermediates)
            summary.update(mt_code=code,best_candidate=candidate,annotation_ok=result['ok'],annotation_status=result['status'],
                           annotation_scope=result['scope'],annotation_candidate_id=result['candidate_id'],
                           annotation_dir=str(sample/'MitoLooper') if result['ok'] else result.get('directory',''))
            summary['status']=recovered_status(summary)
            summary['postprocessed_at']=datetime.now(timezone.utc).isoformat(timespec='seconds')
            changed={event['stage'] for event in runner.events}
            stages=[event for event in old_stages if event.get('stage') not in changed]+runner.events
            summary['failed_stages']=[event['stage'] for event in stages if event.get('status') in {'failed','timeout'}]
            c.atomic_json(sample/'stages.json',stages)
            m.write_table(sample/'stages.tsv',stages,['stage','status','reason','returncode','elapsed_seconds','log','command'])
            c.atomic_json(path,summary)
            candidates_path=sample/'candidates.tsv'
            if candidates_path.is_file():
                import csv
                with candidates_path.open(newline='') as handle:
                    reader=csv.DictReader(handle,delimiter='\t')
                    fields=list(reader.fieldnames or [])
                    candidates=list(reader)
                if 'target_reason' not in fields:
                    fields.append('target_reason')
                for row in candidates:
                    if row.get('candidate_id')==candidate['candidate_id']:
                        row.update({key:value for key,value in candidate.items() if key in fields})
                m.write_table(candidates_path,candidates,fields)
            c.standardize_layout(sample)
            c.atomic_json(cache,{'signature':signature,'status':result['status'],'scope':result['scope']})
            rows.append({'sample':sample.name,'status':result['status'],'bold_recheck_failed':bool(recheck_error)})
        except (ValueError,OSError,RuntimeError) as error:
            print(f'  FAILED: {error}',flush=True)
            rows.append({'sample':sample.name,'status':'failed','reason':str(error)})
    c.refresh_cohort(root)
    c.atomic_json(root/'postprocess_summary.json',{'results':rows,'rechecked_bold':cli.recheck_bold})
    if not cli.keep_intermediates:
        directory=root/'postprocess_bold_reference'
        if directory.is_dir() and not directory.is_symlink():
            for path in c.owned_files(directory):
                if path.name.startswith('coi.'):
                    path.unlink()
    return rows


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('run_dir',type=Path)
    parser.add_argument('--recheck-bold',action='store_true',help='Repeat COI validation using saved localization when available and the saved BOLD inputs')
    parser.add_argument('--only-sample')
    parser.add_argument('--force',action='store_true',help='Repeat an otherwise cached annotation')
    parser.add_argument('--keep-intermediates',action='store_true')
    parser.add_argument('--mt-code',type=int,help='Explicit code override; otherwise use each saved sample code')
    parser.add_argument('--mitos-refdir',type=Path)
    parser.add_argument('--mitos-refseqver')
    parser.add_argument('--bold-fasta',type=Path)
    parser.add_argument('--bold-metadata',type=Path)
    parser.add_argument('--bold-max-hits',type=int)
    parser.add_argument('--bold-search-max-hits',type=int)
    parser.add_argument('--threads',type=int)
    parser.add_argument('--annotation-timeout',type=float)
    parser.add_argument('--config',type=Path)
    cli=parser.parse_args(argv)
    root=cli.run_dir.absolute()
    if root.is_symlink() or not root.is_dir():
        parser.error('Not a regular Tortuga run directory')
    with (root/'.pipeline.lock').open('a') as lock:
        try:
            fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError:
            parser.error('The pipeline is still active; no annotation or validation was started')
        rows=postprocess(root,cli)
    if not rows:
        parser.error('No completed mitochondrial samples with retained sequences were selected')
    return 1 if any((row['status'] not in {'completed','cached','skipped_no_candidate'} or row.get('bold_recheck_failed')) for row in rows) else 0


if __name__=='__main__':
    raise SystemExit(main())
