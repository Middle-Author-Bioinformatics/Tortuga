'use strict';
const report = JSON.parse(document.getElementById('report-data').textContent);
const state = {sample: 0, branch: 'mito', mode: 'map', locus: 0};
const $ = id => document.getElementById(id);
const fmt = value => value === null || value === undefined || value === '' ? '—' : Number(value).toLocaleString(undefined, {maximumFractionDigits: 1});
const human = value => String(value || 'unknown').replaceAll('_', ' ');
const methodName = value => ({novoplasty:'NOVOPlasty',getorganelle:'GetOrganelle',getorganelle_assembly:'GetOrganelle',mitofinder:'MitoFinder',mitofinder_assembly:'MitoFinder',mtgrasp:'mtGrasp'})[value] || human(value);
const exists = value => value !== null && value !== undefined && value !== '';
function el(tag, attrs = {}, ...children) {
  const node = document.createElement(tag);
  for (const [key, value] of Object.entries(attrs)) {
    if (key.startsWith('on')) node.addEventListener(key.slice(2), value);
    else if (value !== null && value !== undefined) node.setAttribute(key, value);
  }
  for (const child of children.flat()) if (child !== null && child !== undefined) node.append(child instanceof Node ? child : document.createTextNode(String(child)));
  return node;
}
function svg(tag, attrs = {}, ...children) {
  const node = document.createElementNS('http://www.w3.org/2000/svg', tag);
  for (const [key, value] of Object.entries(attrs)) node.setAttribute(key, value);
  children.flat().forEach(child => node.append(child));
  return node;
}
const empty = (title, message) => el('div', {class: 'empty'}, el('strong', {}, title), message);
function link(file, label, primary = false) {
  return file ? el('a', {href: file.href, download: file.name, class: primary ? 'primary-link' : ''}, label || file.name) : null;
}
function metrics(items) {
  return el('div', {class: 'metrics'}, items.map(([label, value, suffix, note]) => el('div', {class: 'metric'}, el('span', {}, label), el('strong', {}, value, suffix ? el('small', {}, suffix) : null), note ? el('p', {}, note) : null)));
}
function table(headers, rows) {
  return el('div', {class: 'table-wrap'}, el('table', {}, el('thead', {}, el('tr', {}, headers.map(header => el('th', {scope: 'col'}, header)))), el('tbody', {}, rows.map(row => el('tr', {}, row.map(cell => el('td', {}, exists(cell) ? cell : '—')))))));
}
function details(title, content) {return el('details', {class: 'details'}, el('summary', {}, title), content);}
let sequencePopover;
function showSequence(anchor, name, sequence, description = '') {
  sequencePopover?.remove();
  const field = el('textarea', {readonly: '', 'aria-label': 'Nucleotide sequence', spellcheck: 'false'}, sequence || 'Sequence unavailable');
  const copy = el('button', {class: 'copy-sequence', disabled: sequence ? null : '', onclick: async () => {
    try {await navigator.clipboard.writeText(sequence); copy.textContent = 'Copied';}
    catch (_) {field.focus(); field.select(); if (document.execCommand) document.execCommand('copy'); copy.textContent = 'Selected';}
  }}, 'Copy sequence');
  const fasta = '>'+name.replaceAll(/[^A-Za-z0-9_.-]/g, '_')+'\n'+(sequence.match(/.{1,80}/g) || []).join('\n')+'\n';
  const box = el('div', {class: 'sequence-popover', role: 'dialog', 'aria-label': name+' sequence'},
    el('div', {class: 'popover-head'}, el('strong', {}, name), el('button', {'aria-label': 'Close sequence', onclick: () => box.remove()}, '×')),
    el('p', {}, description), field,
    el('div', {class: 'file-line'}, copy, sequence ? el('a', {href: 'data:text/plain;charset=utf-8,'+encodeURIComponent(fasta), download: name+'.fasta'}, 'FASTA ↓') : null,
      el('a', {href: 'https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch', target: '_blank', rel: 'noopener noreferrer'}, 'NCBI BLAST ↗')));
  document.body.append(box); sequencePopover = box;
  const bounds = anchor.getBoundingClientRect();
  box.style.left = Math.max(12, Math.min(bounds.left, window.innerWidth-440))+'px';
  box.style.top = Math.max(12, Math.min(bounds.bottom+8, window.innerHeight-350))+'px';
}
function fileList(files, title = 'All retained files') {
  return files.length ? details(title+' · '+files.length, table(['File', 'Size'], files.map(file => [link(file), fmt(file.bytes/1024)+' KiB']))) : null;
}
function selectSample(index) {state.sample = index; state.locus = 0; state.mode = 'map'; render();}
function sampleList() {
  const query = $('search').value.toLowerCase();
  const list = $('sample-list'); list.replaceChildren();
  report.samples.forEach((sample, index) => {
    const text = [sample.id, sample.expected_taxon, ...(sample.mito.reference_taxa || []), ...sample.its.loci.map(locus => locus.taxon)].join(' ').toLowerCase();
    if (!text.includes(query)) return;
    const mitoClass = sample.mito.accepted ? 'yes' : sample.mito.sequence ? 'partial' : sample.status === 'failed' ? 'no' : '';
    const itsClass = sample.its.count > 0 ? 'yes' : '';
    list.append(el('button', {class: 'sample-row', 'aria-pressed': String(index === state.sample), onclick: () => selectSample(index)},
      el('span', {}, el('strong', {}, sample.id), el('small', {}, sample.expected_taxon || sample.its.loci[0]?.taxon || 'Identity not supplied')),
      el('span', {class: 'dots', 'aria-label': `Mitogenome ${sample.mito.accepted ? 'accepted' : 'not accepted'}; ${sample.its.count || 0} complete ITS loci`},
        el('span', {class: 'dot '+mitoClass, 'aria-hidden': 'true'}), el('span', {class: 'dot '+itsClass, 'aria-hidden': 'true'}))));
  });
  if (!list.children.length) list.append(el('p', {class: 'blank-list'}, 'No matching samples.'));
}
function viewSwitch(labels) {
  return el('div', {class: 'view-switch', 'aria-label': 'Assembly display'}, labels.map(([key, label]) => el('button', {class: state.mode === key ? 'active' : '', 'aria-pressed': String(state.mode === key), onclick: () => {state.mode = key; renderView();}}, label)));
}
function mapSequence(length, features, regionMap = false) {
  if (!length || !features.length) return empty(regionMap ? 'No region coordinates available' : 'No annotation map available', '');
  const panel = el('div', {class: 'sequence-panel'});
  const drawing = svg('svg', {class: 'genome-map', viewBox: '0 0 940 174', role: 'group', 'aria-label': regionMap ? 'ITS locus regions on a base-pair coordinate axis' : 'Mitochondrial features on an unrolled base-pair coordinate axis'});
  const scale = position => 45 + (position-1)/length*850;
  drawing.append(svg('line', {x1: 45, x2: 895, y1: 88, y2: 88, stroke: '#486053', 'stroke-width': 1}));
  for (let i = 0; i <= 4; i++) {
    const position = Math.round(i*length/4), x = 45+i*850/4;
    drawing.append(svg('line', {x1: x, x2: x, y1: 27, y2: 145, stroke: '#253c32', 'stroke-width': 1, 'stroke-dasharray': '2 5'}));
    const tick = svg('text', {x, y: 164, class: 'tick', 'text-anchor': 'middle'}); tick.textContent = i === 0 ? '1' : fmt(position)+' bp'; drawing.append(tick);
  }
  if (!regionMap) {
    for (const [label, y] of [['+', 67], ['−', 119]]) {const text = svg('text', {x: 20, y, class: 'tick'}); text.textContent = label; drawing.append(text);}
  }
  const readout = el('div', {class: 'feature-readout', 'aria-live': 'polite'}, 'Select a feature');
  const colors = regionMap ? {'SSU': '#587c69', 'ITS1': '#96d1b1', '5.8S': '#d4b477', 'ITS2': '#96d1b1', 'LSU': '#587c69'} : {CDS: '#96d1b1', rRNA: '#d4b477', tRNA: '#b7a1ce'};
  for (const feature of features) {
    const name = regionMap ? feature.name : feature.gene;
    const y = regionMap ? 65 : feature.strand === '-' ? 99 : 47;
    const group = svg('g', {class: 'gene', tabindex: 0, role: 'button', 'aria-label': `${name}, ${feature.start} to ${feature.end} bp`});
    for (const [start, end] of feature.segments || [[feature.start, feature.end]]) {
      const x = scale(start), width = Math.max(2, (end-start+1)/length*850);
      group.append(svg('rect', {x, y, width, height: 29, rx: 2, fill: colors[regionMap ? name : feature.type] || '#7e9a8a'}));
      if (width > name.length*6.7+8) {const label = svg('text', {x: x+width/2, y: y+19, 'text-anchor': 'middle', fill: '#10251b', 'font-size': 11}); label.textContent = name; group.append(label);}
    }
    const select = () => {drawing.querySelectorAll('.selected').forEach(node => node.classList.remove('selected')); group.classList.add('selected'); readout.textContent = `${name} · ${fmt(feature.start)}–${fmt(feature.end)} bp${feature.strand ? ' · '+feature.strand : ''}`; showSequence(group, name, feature.sequence || '', readout.textContent);};
    group.addEventListener('click', select); group.addEventListener('keydown', event => {if (event.key === 'Enter' || event.key === ' ') {event.preventDefault(); select();}});
    drawing.append(group);
  }
  panel.append(drawing, readout, el('div', {class: 'legend'}, Object.entries(regionMap ? {'ITS cassette': '#96d1b1', '5.8S': '#d4b477', 'Flanking rRNA': '#587c69'} : {CDS: '#96d1b1', rRNA: '#d4b477', tRNA: '#b7a1ce'}).map(([label, color]) => el('span', {}, el('i', {style: 'background:'+color}), label))));
  return panel;
}
function renderGraph(graphs, host, locus = null, method = '') {
  if (locus) graphs = graphs.filter(graph => (graph.loci || []).includes(locus.locus));
  if (!graphs.length) {host.append(empty('Source graph unavailable', method === 'novoplasty' ? 'NOVOPlasty sequence output' : '')); return;}
  const select = el('select', {'aria-label': 'Assembly graph file'}, graphs.map((graph, index) => el('option', {value: index}, graph.label || graph.file?.name || 'Selected assembly')));
  const zoomIn = el('button', {'aria-label': 'Zoom in'}, '+'), zoomOut = el('button', {'aria-label': 'Zoom out'}, '−'), reset = el('button', {'aria-label': 'Fit graph'}, 'Fit');
  const tools = el('div', {class: 'graph-tools'}, select, el('div', {class: 'graph-buttons'}, zoomOut, zoomIn, reset));
  const plot = el('div'); host.append(tools, plot);
  const draw = () => {
    const graph = graphs[Number(select.value)]; plot.replaceChildren();
    [zoomIn, zoomOut, reset].forEach(button => {button.disabled = !graph.nodes.length;});
    if (!graph.nodes.length) {plot.append(empty('Graph preview unavailable', ''), el('div', {class: 'file-line'}, link(graph.file, 'Download graph'))); return;}
    const drawing = svg('svg', {class: 'graph', viewBox: '0 0 1000 370', role: 'group', 'aria-label': 'Interactive assembly graph. Scroll to zoom, drag to pan, select a segment for details.'});
    const pane = svg('g'); drawing.append(pane);
    const lookup = new Map(graph.nodes.map(node => [node.id, node]));
    for (const edge of graph.edges) {
      const a = lookup.get(edge.source), b = lookup.get(edge.target);
      const start = a.x+(edge.source_strand === '+' ? 34 : -34), end = b.x+(edge.target_strand === '+' ? -34 : 34);
      const path = edge.source === edge.target ? `M ${start} ${a.y} C ${a.x+100} ${a.y-100}, ${a.x-100} ${a.y-100}, ${end} ${b.y}` : `M ${start} ${a.y} C ${(start+end)/2} ${a.y}, ${(start+end)/2} ${b.y}, ${end} ${b.y}`;
      const curve = svg('path', {class: 'link', d: path}); const title = svg('title'); title.textContent = `${edge.source}${edge.source_strand} → ${edge.target}${edge.target_strand}${edge.overlap ? ' · overlap '+edge.overlap : ''}`; curve.append(title); pane.append(curve);
    }
    const readout = el('div', {class: 'graph-readout', 'aria-live': 'polite', title: 'Graph positions and segment widths are schematic.'}, 'Drag to pan · scroll to zoom · select a segment');
    const selectedLocus = locus?.locus || 'mitogenome';
    for (const node of graph.nodes) {
      const active = (node.loci || []).includes(selectedLocus);
      const group = svg('g', {class: 'segment'+(active?' path-node':' context-node')+((node.loci || []).length > 1?' shared-node':''), transform: `translate(${node.x} ${node.y})`, tabindex: 0, role: 'button', 'aria-label': 'Segment '+node.id});
      group.append(svg('rect', {x: -34, y: -13, width: 68, height: 26, rx: 8}));
      const text = svg('text', {'text-anchor': 'middle', y: 4}); text.textContent = node.id.length > 9 ? node.id.slice(0, 8)+'…' : node.id; group.append(text);
      const taxa = node.taxa || [];
      if (taxa.length) {const label = svg('text', {'text-anchor': 'middle', y: 30, class: 'node-taxon'}); label.textContent = taxa.length > 1 ? taxa.length+' taxa' : taxa[0].length > 20 ? taxa[0].slice(0, 19)+'…' : taxa[0]; group.append(label);}
      const choose = () => {pane.querySelectorAll('.chosen').forEach(item => item.classList.remove('chosen')); group.classList.add('chosen'); readout.textContent = `Segment ${node.id} · ${fmt(node.length)} bp · ${exists(node.coverage) ? fmt(node.coverage)+'×' : '—'} · ${taxa.join(' / ') || 'Taxonomy unavailable'}`; showSequence(group, 'Segment_'+node.id, node.sequence || '', `${readout.textContent} · ${node.taxonomy_scope || ''} · ${(node.loci || []).join(' / ')}`);};
      group.addEventListener('click', choose); group.addEventListener('keydown', event => {if (event.key === 'Enter' || event.key === ' ') {event.preventDefault(); choose();}}); pane.append(group);
    }
    let tx = 0, ty = 0, scale = 1;
    const transform = () => pane.setAttribute('transform', `translate(${tx} ${ty}) scale(${scale})`);
    const fit = () => {const xs = graph.nodes.map(node => node.x), ys = graph.nodes.map(node => node.y); const minX = Math.min(...xs)-100, maxX = Math.max(...xs)+100, minY = Math.min(...ys)-100, maxY = Math.max(...ys)+55; scale = Math.min(900/(maxX-minX), 290/(maxY-minY), 1.8); tx = (1000-(maxX-minX)*scale)/2-minX*scale; ty = (370-(maxY-minY)*scale)/2-minY*scale; transform();};
    const zoom = (factor, x = 500, y = 185) => {const next = Math.max(.001, Math.min(12, scale*factor)); const ratio = next/scale; tx = x-(x-tx)*ratio; ty = y-(y-ty)*ratio; scale = next; transform();};
    zoomIn.onclick = () => zoom(1.35); zoomOut.onclick = () => zoom(1/1.35); reset.onclick = fit;
    const point = event => {const bounds = drawing.getBoundingClientRect(); const matrix = drawing.getScreenCTM(); return matrix ? new DOMPoint(event.clientX, event.clientY).matrixTransform(matrix.inverse()) : {x: (event.clientX-bounds.left)*1000/bounds.width, y: (event.clientY-bounds.top)*370/bounds.height};};
    drawing.addEventListener('wheel', event => {event.preventDefault(); const p = point(event); zoom(event.deltaY < 0 ? 1.15 : 1/1.15, p.x, p.y);}, {passive: false});
    let drag;
    drawing.addEventListener('pointerdown', event => {if (event.target.closest('.segment')) return; drag = point(event); drawing.setPointerCapture(event.pointerId);});
    drawing.addEventListener('pointermove', event => {if (!drag) return; const p = point(event); tx += p.x-drag.x; ty += p.y-drag.y; drag = p; transform();});
    drawing.addEventListener('pointerup', () => {drag = null;}); drawing.addEventListener('pointercancel', () => {drag = null;});
    fit();
    plot.append(drawing, readout, el('div', {class: 'graph-caption'}, el('span', {}, `${fmt(graph.nodes.length)}${graph.truncated ? ' / '+fmt(graph.total_nodes) : ''} segments · ${fmt(graph.edges.length)} links`), link(graph.file, 'Download graph')));
    const path = graph.paths?.[selectedLocus];
    if (path) {
      const sum = path.steps.reduce((total, step) => total+step.length, 0);
      const overlaps = path.steps.reduce((total, step) => total+(step.overlap || 0), 0);
      plot.append(el('div', {class: 'path-ribbon'}, path.steps.map(step => el('span', {}, step.node+step.orientation))));
      if (path.sequence_verified) plot.append(el('p', {class: 'path-equation'}, `${fmt(path.locus_bp)} bp = ${fmt(sum)} segment bp − ${fmt(overlaps)} overlap bp${path.trimmed_bp ? ' − '+fmt(path.trimmed_bp)+' trimmed bp' : ''}`));
      plot.append(details('Sequence path', table(['Segment', 'Strand', 'Length', 'Overlap', 'Added bp'], path.steps.map(step => [step.node, step.orientation, fmt(step.length), fmt(step.overlap), fmt(step.added_bp)]))));
    }
  };
  select.addEventListener('change', draw); draw();
}
function mitoView(sample) {
  const mito = sample.mito, root = document.createDocumentFragment();
  root.append(metrics([
    [mito.accepted ? 'Accepted assembly' : 'Best candidate', fmt(mito.length), 'bp'],
    ['Topology', mito.circular === true ? 'Circular' : mito.circular === false ? 'Unclosed' : 'Unknown'],
    [mito.core_genes?.source === 'MITOS annotation' ? 'Annotated core CDS / rRNA' : 'Screened core CDS / rRNA', `${fmt(mito.core_genes?.cds)} / 13 · ${fmt(mito.core_genes?.rrna)} / 2`, '', mito.core_genes?.extra?.length ? 'Additional: '+mito.core_genes.extra.join(', ') : ''],
    ['Recovery method', mito.method ? methodName(mito.method) : '—']
  ]));
  root.append(el('div', {class: 'section-head'}, el('h3', {}, state.mode === 'graph' ? 'Assembly topology' : 'Mitochondrial sequence'), viewSwitch([['map', 'Sequence & genes'], ['graph', 'Assembly graph']])));
  const host = el('div'); root.append(host);
  if (state.mode === 'graph') renderGraph(mito.graphs, host, null, mito.method);
  else {
    const features = mito.annotation.features;
    const types = new Set(['CDS', 'gene', 'rRNA', 'tRNA', 'ncRNA']);
    const selected = features.some(feature => types.has(feature.type)) ? features.filter(feature => types.has(feature.type)) : features.filter(feature => feature.type === 'gene');
    host.append(mito.annotation.multiple_sequences ? empty('Multiple annotation sequences', '') : mapSequence(mito.length, selected));
  }
  const supported = mito.target_status === 'confirmed_target';
  const taxa = mito.reference_taxa || [];
  const bestHit = mito.hits[0];
  const confidence = mito.confidence || {};
  root.append(el('div', {class: 'evidence'}, el('div', {}, el('h4', {}, confidence.taxon && !taxa.includes(confidence.taxon) ? 'Shared COI reference lineage' : 'COI reference match'), el('em', {}, confidence.taxon || (taxa.length ? taxa.slice(0, 3).join(' / ')+(taxa.length > 3 ? ' / …' : '') : mito.hits.length ? 'Taxonomy unavailable' : 'No COI match')), bestHit ? el('small', {}, `${fmt(bestHit.identity)}% identity · ${fmt(bestHit.aligned_bp)} bp${exists(bestHit.coverage)?' · '+fmt(bestHit.coverage*100)+'% coverage':''}`) : null), el('div', {}, el('h4', {}, 'Match confidence'), el('p', {title: 'Evidence grade from identity, aligned bases, coverage, metadata and near-best taxonomic agreement; not a probability.'}, confidence.grade || 'Unavailable'))));
  root.append(el('div', {class: 'file-line'}, link(mito.sequence, mito.accepted ? 'Download accepted FASTA ↓' : 'Download candidate FASTA ↓', true), link(mito.annotation.file, 'MITOS annotations'), link(sample.summary_file, 'Evidence summary')));
  if (mito.annotation.features.length) root.append(details(`Annotation table · ${mito.annotation.features.length} feature records`, table(['Gene / feature', 'Type', 'Start', 'End', 'Strand', 'Product'], mito.annotation.features.map(feature => [feature.gene, feature.type, fmt(feature.start), fmt(feature.end), feature.strand, feature.product]))));
  if (mito.reference_genes?.length) root.append(details('Reference gene recovery', table(['Gene', 'Type', 'Detected', 'Coverage'], mito.reference_genes.map(gene => [gene.gene, gene.type, String(gene.screen_pass).toLowerCase() === 'true' ? 'Yes' : 'No', fmt(Number(gene.query_coverage)*100)+'%']))));
  root.append(details('Taxonomic evidence', el('div', {},
    table(['Check', 'Result'], [['Target agreement', !mito.hits.length ? 'Unavailable' : supported ? 'Confirmed' : mito.target_status === 'non_target' ? 'Different taxon' : 'Unconfirmed']]),
    mito.hits.length ? table(['BOLD record', 'Reference taxon', 'Identity', 'Aligned bp'], mito.hits.map(hit => [hit.record, hit.taxon, fmt(hit.identity)+'%', fmt(hit.aligned_bp)])) : null,
    mito.hit_count > mito.hits.length ? el('p', {class: 'subtle'}, `${mito.hits.length} of ${fmt(mito.hit_count)} records shown`) : null,
    el('div', {class: 'file-line'}, link(mito.raw_hits_file, 'COI alignments ↓'), link(mito.evidence_file, 'BOLD results ↓')))));
  if (mito.files.length) root.append(fileList(mito.files));
  return root;
}
function itsView(sample) {
  const its = sample.its, root = document.createDocumentFragment();
  const index = Math.min(state.locus, Math.max(0, its.loci.length-1)), locus = its.loci[index];
  root.append(metrics([
    ['Complete loci', fmt(its.count)],
    ['ITS cassette', locus ? fmt(locus.its_length_bp) : '—', 'bp'],
    ['Full rDNA locus', locus ? fmt(locus.locus_length_bp) : '—', 'bp'],
    ['Consensus resolution', locus?.taxonomy_rank ? human(locus.taxonomy_rank) : 'Unknown']
  ]));
  if (its.loci.length) {
    const select = el('select', {'aria-label': 'Reported ITS locus', onchange: event => {state.locus = Number(event.target.value); renderView();}}, its.loci.map((entry, i) => el('option', {value: i, selected: i === index ? '' : null}, entry.locus+(entry.taxon?' · '+entry.taxon:'')+' · '+fmt(entry.locus_length_bp)+' bp')));
    root.append(el('label', {class: 'locus-select'}, 'Reported locus', select));
  }
  root.append(el('div', {class: 'section-head'}, el('h3', {}, state.mode === 'graph' ? 'rDNA assembly graph' : 'Reconstructed rDNA locus'), viewSwitch([['map', 'Locus map'], ['graph', 'Assembly graph']])));
  const host = el('div'); root.append(host);
  if (state.mode === 'graph') renderGraph(its.graphs, host, locus);
  else if (locus) host.append(mapSequence(locus.locus_length_bp || Math.max(0, ...locus.regions.map(region => region.end)), locus.regions, true));
  else host.append(empty(its.status === 'not_requested' ? 'ITSME not requested' : ['missing_report','unsupported_report','failed_no_report','not_run'].includes(its.status) ? 'ITSME results unavailable' : 'No complete ITS locus reported', ''));
  if (locus) {
    const selection = locus.selection || {};
    const regionNames = new Set(locus.regions.map(region => region.name));
    root.append(details('Locus selection', table(['Criterion', 'Result'], [
      ['Structure', ['SSU','ITS1','5.8S','ITS2','LSU'].filter(name => regionNames.has(name)).join(' → ') || 'Not recorded'],
      ['Taxonomy', String(locus.taxonomy_status).includes('PASS') ? 'Concordant' : 'Tentative'],
      ['Source', selection.source || locus.locus],
      ['Junctions supported', Number.isFinite(Number(selection.supported_junctions)) ? selection.supported_junctions+' / '+selection.junction_count : 'Not recorded'],
      ['Representative', selection.alternative_count ? 'Longest of '+(selection.alternative_count+1)+' validated paths' : 'Reported complete locus']
    ])));
    root.append(el('div', {class: 'evidence'}, el('div', {}, el('h4', {}, 'Predicted consensus taxon'), el('em', {}, locus.taxon || 'Not assigned')), el('div', {}, el('h4', {}, 'Cassette lengths'), el('p', {}, `ITS1 ${fmt(locus.its1_length_bp)} bp · 5.8S ${fmt(locus.five_8s_length_bp)} bp · ITS2 ${fmt(locus.its2_length_bp)} bp`))));
    if (locus.consensus_taxonomy) root.append(el('p', {class: 'lineage'}, el('strong', {}, 'Consensus lineage  '), locus.consensus_taxonomy));
  }
  root.append(el('div', {class: 'file-line'}, link(its.sequence, its.sequence?.name.includes('complete_ITS') ? 'ITS cassette FASTA ↓' : 'rDNA loci FASTA ↓', true), its.cassette_sequence?.href !== its.sequence?.href ? link(its.cassette_sequence, 'ITS cassette FASTA') : null, link(its.summary_file, 'ITSME results table')));
  if (its.loci.length) root.append(details('All reported complete loci · '+its.loci.length, table(['Locus', 'ITS bp', 'Full locus bp', 'Consensus taxon', 'Rank'], its.loci.map(entry => [entry.locus, fmt(entry.its_length_bp), fmt(entry.locus_length_bp), entry.taxon, entry.taxonomy_rank]))));
  if (its.files.length) root.append(fileList(its.files)); return root;
}
function runView(sample) {
  const root = document.createDocumentFragment(), qc = sample.qc;
  root.append(metrics([['Input reads', fmt(qc.reads_before)], ['Retained reads', fmt(qc.reads_after), '', exists(qc.reads_before) && qc.reads_before > 0 && exists(qc.reads_after) ? fmt(qc.reads_after/qc.reads_before*100)+'% retained' : ''], ['Q30 after QC', exists(qc.q30_after) ? fmt(qc.q30_after*100)+'%' : '—'], ['Recorded stages', fmt(sample.stages.length)]]));
  root.append(el('div', {class: 'progress-track', 'aria-hidden': 'true'}, sample.stages.map(stage => el('div', {class: stage.status.toLowerCase(), title: stage.stage+': '+stage.status}))));
  if (sample.stages.length) root.append(table(['Stage', 'Status', 'Time'], sample.stages.map(stage => {
    const failed = ['failed', 'timeout'].includes(stage.status.toLowerCase());
    const status = failed && stage.log ? link(stage.log, human(stage.status)+' ↓') : el('span', {}, human(stage.status));
    status.className = 'stage-status '+stage.status.toLowerCase();
    if (failed && stage.log) status.setAttribute('aria-label', `Download ${human(stage.stage)} ${human(stage.status)} log`);
    return [human(stage.stage), status, numberTime(stage.elapsed_seconds)];
  })));
  else root.append(empty('No stage records', ''));
  root.append(el('div', {class: 'file-line'}, link(qc.report, 'Open QC report'), link(sample.summary_file, 'Sample summary')));
  if (sample.logs.length) root.append(fileList(sample.logs, 'Native logs'));
  return root;
}
function numberTime(value) {const n = Number(value); return !exists(value) || !Number.isFinite(n) ? '—' : n < 60 ? fmt(n)+' s' : fmt(n/60)+' min';}
function renderView() {
  sequencePopover?.remove();
  const sample = report.samples[state.sample];
  $('view').replaceChildren(state.branch === 'mito' ? mitoView(sample) : state.branch === 'its' ? itsView(sample) : runView(sample));
  const badge=$('sample-status');
  badge.textContent=state.branch==='its'?'ITSME · '+human(sample.its.run_status):state.branch==='mito'?'Mito · '+(sample.mito.accepted?'Accepted':sample.mito.sequence?'Candidate':'No sequence'):'Run record';
  badge.className='status '+(state.branch==='its'?(sample.its.run_status==='completed'?'':['failed','timeout'].includes(sample.its.run_status)?'fail':'warn'):sample.status.includes('failed')?'fail':sample.mito.accepted?'':'warn');
  document.querySelectorAll('[data-branch]').forEach(button => {if (button.dataset.branch === state.branch) button.setAttribute('aria-current', 'page'); else button.removeAttribute('aria-current');});
}
function render() {
  const sample = report.samples[state.sample]; sampleList();
  $('sample-name').textContent = sample.id;
  $('sample-position').textContent = `SAMPLE ${String(state.sample+1).padStart(2, '0')} / ${String(report.samples.length).padStart(2, '0')}`;
  $('expected-name').replaceChildren(sample.expected_taxon ? 'Expected identity · ' : 'Expected identity not supplied', sample.expected_taxon ? el('em', {}, sample.expected_taxon) : '');
  const ranks = ['domain','kingdom','phylum','class','order','family','genus','species'];
  const hierarchy = Object.entries(sample.expected_taxonomy || {}).filter(([rank, name]) => name).sort(([a],[b]) => ranks.indexOf(a.toLowerCase())-ranks.indexOf(b.toLowerCase())).map(([rank, name]) => el('span', {title: human(rank)}, name));
  $('expected-hierarchy').replaceChildren(...hierarchy.flatMap((item, i) => i ? [' › ', item] : [item]));
  const badge = $('sample-status'); badge.textContent = human(sample.status); badge.className = 'status '+(sample.status.includes('failed') ? 'fail' : sample.mito.accepted ? '' : 'warn');
  renderView();
}
function compare() {
  const host = $('comparison-chart'); host.replaceChildren();
  const maxMito = Math.max(1, ...report.samples.map(sample => sample.mito.length || 0));
  const maxIts = Math.max(1, ...report.samples.flatMap(sample => sample.its.loci.map(locus => locus.its_length_bp || 0)));
  host.append(el('div', {class: 'compare-row compare-head'}, el('span', {}, 'Sample'), el('span', {}, 'MITOCHONDRIAL BEST · BP'), el('span', {}, 'LONGEST COMPLETE ITS · BP')));
  report.samples.forEach((sample, index) => {
    const mito = sample.mito.length, lengths = sample.its.loci.map(locus => locus.its_length_bp).filter(exists), its = lengths.length ? Math.max(...lengths) : null;
    host.append(el('div', {class: 'compare-row', tabindex: 0, role: 'button', 'aria-label': 'Open sample '+sample.id,
      onclick: () => {$('comparison').close(); selectSample(index);}, onkeydown: event => {if (event.key === 'Enter' || event.key === ' ') {event.preventDefault(); $('comparison').close(); selectSample(index);}}},
      el('span', {}, sample.id), el('div', {}, fmt(mito), el('small', {}, sample.mito.accepted ? ' · accepted' : ' · candidate'), mito ? el('div', {class: 'compare-bar', style: `width:${mito/maxMito*100}%`}) : null),
      el('div', {}, fmt(its), its ? el('div', {class: 'compare-bar its', style: `width:${its/maxIts*100}%`}) : null)));
  });
  $('comparison').showModal();
}
$('report-title').textContent = report.title;
document.title = report.title+' · Tortuga';
$('cohort-counts').textContent = `${report.samples.length} samples · ${report.samples.filter(sample => sample.mito.accepted).length} mito accepted · ${report.samples.filter(sample => sample.its.count > 0).length} with ITS`;
$('built').textContent = 'Built '+report.built.slice(0, 10)+' · UTC';
$('cohort-files').append(...report.files.slice(0, 2).map(file => link(file, file.name === 'cohort_summary.csv' ? 'Cohort CSV ↓' : 'ITS loci CSV ↓')));
$('search').addEventListener('input', sampleList);
$('branch-tabs').addEventListener('click', event => {const button = event.target.closest('[data-branch]'); if (button) {state.branch = button.dataset.branch; state.mode = 'map'; renderView();}});
$('compare').addEventListener('click', compare); $('close-comparison').addEventListener('click', () => $('comparison').close());
render();
