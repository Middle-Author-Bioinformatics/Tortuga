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
function fileList(files, title = 'All retained files') {
  return files.length ? details(title+' · '+files.length, table(['File', 'Size'], files.map(file => [link(file), fmt(file.bytes/1024)+' KiB']))) : null;
}
function selectSample(index) {state.sample = index; state.locus = 0; state.mode = 'map'; render();}
function sampleList() {
  const query = $('search').value.toLowerCase();
  const list = $('sample-list'); list.replaceChildren();
  report.samples.forEach((sample, index) => {
    const text = [sample.id, sample.expected_taxon, ...sample.its.loci.map(locus => locus.taxon)].join(' ').toLowerCase();
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
  if (!length || !features.length) return empty(regionMap ? 'No region coordinates available' : 'No annotation map available', regionMap ? 'Use the retained ITSME sequences and taxonomy tables below.' : 'The retained FASTA may still be available. MITOS annotations are shown when produced.');
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
  const readout = el('div', {class: 'feature-readout', 'aria-live': 'polite'}, 'Select a feature to see its coordinates.');
  const colors = regionMap ? {'SSU': '#587c69', 'ITS1': '#96d1b1', '5.8S': '#d4b477', 'ITS2': '#96d1b1', 'LSU': '#587c69'} : {CDS: '#96d1b1', rRNA: '#d4b477', tRNA: '#b7a1ce'};
  for (const feature of features) {
    const name = regionMap ? feature.name : feature.gene;
    const x = scale(feature.start), width = Math.max(2, (feature.end-feature.start+1)/length*850);
    const y = regionMap ? 65 : feature.strand === '-' ? 99 : 47;
    const group = svg('g', {class: 'gene', tabindex: 0, role: 'button', 'aria-label': `${name}, ${feature.start} to ${feature.end} bp`});
    group.append(svg('rect', {x, y, width, height: 29, rx: 2, fill: colors[regionMap ? name : feature.type] || '#7e9a8a'}));
    if (width > name.length*6.7+8) {const label = svg('text', {x: x+width/2, y: y+19, 'text-anchor': 'middle', fill: '#10251b', 'font-size': 11}); label.textContent = name; group.append(label);}
    const select = () => {drawing.querySelectorAll('.selected').forEach(node => node.classList.remove('selected')); group.classList.add('selected'); readout.textContent = `${name} · ${regionMap ? 'region' : feature.type} · ${fmt(feature.start)}–${fmt(feature.end)} bp · ${fmt(feature.end-feature.start+1)} bp${feature.strand ? ' · strand '+feature.strand : ''}${feature.product ? ' · '+feature.product : ''}`;};
    group.addEventListener('click', select); group.addEventListener('keydown', event => {if (event.key === 'Enter' || event.key === ' ') {event.preventDefault(); select();}});
    drawing.append(group);
  }
  panel.append(drawing, readout, el('div', {class: 'legend'}, Object.entries(regionMap ? {'ITS cassette': '#96d1b1', '5.8S': '#d4b477', 'Flanking rRNA': '#587c69'} : {CDS: '#96d1b1', rRNA: '#d4b477', tRNA: '#b7a1ce'}).map(([label, color]) => el('span', {}, el('i', {style: 'background:'+color}), label))));
  return panel;
}
function renderGraph(graphs, host) {
  if (!graphs.length) {host.append(empty('No assembly graph retained', 'A final FASTA does not imply that an assembly graph was produced. Available sequence downloads are listed below.')); return;}
  const select = el('select', {'aria-label': 'Assembly graph file'}, graphs.map((graph, index) => el('option', {value: index}, graph.file.name)));
  const zoomIn = el('button', {'aria-label': 'Zoom in'}, '+'), zoomOut = el('button', {'aria-label': 'Zoom out'}, '−'), reset = el('button', {'aria-label': 'Fit graph'}, 'Fit');
  const tools = el('div', {class: 'graph-tools'}, select, el('div', {class: 'graph-buttons'}, zoomOut, zoomIn, reset));
  const plot = el('div'); host.append(tools, plot);
  const draw = () => {
    const graph = graphs[Number(select.value)]; plot.replaceChildren();
    [zoomIn, zoomOut, reset].forEach(button => {button.disabled = !graph.nodes.length;});
    plot.append(el('p', {class: 'subtle'}, graph.scope));
    if (!graph.nodes.length) {plot.append(empty('Graph preview unavailable', graph.warnings.join(' ')), el('div', {class: 'file-line'}, link(graph.file, 'Download original graph'))); return;}
    const drawing = svg('svg', {class: 'graph', viewBox: '0 0 1000 370', role: 'group', 'aria-label': 'Interactive assembly graph. Scroll to zoom, drag to pan, select a segment for details.'});
    const pane = svg('g'); drawing.append(pane);
    const lookup = new Map(graph.nodes.map(node => [node.id, node]));
    for (const edge of graph.edges) {
      const a = lookup.get(edge.source), b = lookup.get(edge.target);
      const start = a.x+(edge.source_strand === '+' ? 34 : -34), end = b.x+(edge.target_strand === '+' ? -34 : 34);
      const path = edge.source === edge.target ? `M ${start} ${a.y} C ${a.x+100} ${a.y-100}, ${a.x-100} ${a.y-100}, ${end} ${b.y}` : `M ${start} ${a.y} C ${(start+end)/2} ${a.y}, ${(start+end)/2} ${b.y}, ${end} ${b.y}`;
      const curve = svg('path', {class: 'link', d: path}); const title = svg('title'); title.textContent = `${edge.source}${edge.source_strand} → ${edge.target}${edge.target_strand}${edge.overlap ? ' · overlap '+edge.overlap : ''}`; curve.append(title); pane.append(curve);
    }
    const readout = el('div', {class: 'graph-readout', 'aria-live': 'polite'}, 'Drag to pan · scroll to zoom · select a segment. Graph positions are a layout, not genomic coordinates.');
    for (const node of graph.nodes) {
      const group = svg('g', {class: 'segment', transform: `translate(${node.x} ${node.y})`, tabindex: 0, role: 'button', 'aria-label': 'Segment '+node.id});
      group.append(svg('rect', {x: -34, y: -13, width: 68, height: 26, rx: 8}));
      const text = svg('text', {'text-anchor': 'middle', y: 4}); text.textContent = node.id.length > 9 ? node.id.slice(0, 8)+'…' : node.id; group.append(text);
      const choose = () => {pane.querySelectorAll('.chosen').forEach(item => item.classList.remove('chosen')); group.classList.add('chosen'); readout.textContent = `Segment ${node.id} · ${fmt(node.length)} bp · coverage ${exists(node.coverage) ? fmt(node.coverage)+'×' : 'not recorded'}. Segment width is schematic; use these values for length.`;};
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
    plot.append(drawing, readout, el('div', {class: 'graph-caption'}, el('span', {}, `${fmt(graph.nodes.length)} displayed segments · ${fmt(graph.edges.length)} oriented links${graph.truncated ? ' · bounded preview' : ''}`), link(graph.file, 'Original graph ↗')));
    graph.warnings.forEach(message => plot.append(el('p', {class: 'note warning'}, message)));
  };
  select.addEventListener('change', draw); draw();
}
function mitoView(sample) {
  const mito = sample.mito, root = document.createDocumentFragment();
  root.append(metrics([
    [mito.accepted ? 'Accepted assembly' : 'Best candidate', fmt(mito.length), 'bp', mito.accepted ? 'Meets the configured stopping policy' : 'May be partial or unresolved'],
    ['Topology', mito.circular === true ? 'Circular' : mito.circular === false ? 'Unclosed' : 'Unknown', '', 'Assembly evidence'],
    ['Reference genes', exists(mito.genes_expected) ? `${fmt(mito.genes_recovered)} / ${fmt(mito.genes_expected)}` : 'Unknown', '', mito.complete_by_reference ? 'Reference-gene screen passed' : 'Homology recovery; not tRNA completeness'],
    ['Recovery method', mito.method ? methodName(mito.method) : '—', '', human(sample.status)]
  ]));
  root.append(el('div', {class: 'section-head'}, el('div', {}, el('h3', {}, state.mode === 'graph' ? 'Assembly topology' : 'The mitochondrial sequence'), el('p', {}, state.mode === 'graph' ? 'Inspect the retained assembly-method graph.' : mito.circular ? 'Circular assembly, shown unrolled from its reported origin.' : 'Reported sequence coordinates.')), viewSwitch([['map', 'Sequence & genes'], ['graph', 'Assembly graph']])));
  const host = el('div'); root.append(host);
  if (state.mode === 'graph') renderGraph(mito.graphs, host);
  else {
    const features = mito.annotation.features;
    const types = new Set(['CDS', 'rRNA', 'tRNA', 'ncRNA']);
    const selected = features.some(feature => types.has(feature.type)) ? features.filter(feature => types.has(feature.type)) : features.filter(feature => feature.type === 'gene');
    host.append(mito.annotation.multiple_sequences ? empty('Multiple annotation sequences', 'See the annotation table below; features are not overlaid across different contigs.') : mapSequence(mito.length, selected));
    if (mito.annotation.file && mito.annotation_ok === false) host.append(el('p', {class: 'note warning'}, 'The annotation stage did not complete successfully. These are retained annotation features; inspect the run record.'));
    mito.annotation.warnings.forEach(message => host.append(el('p', {class: 'note warning'}, message)));
  }
  const supported = mito.target_status === 'confirmed_target';
  root.append(el('div', {class: 'evidence'}, el('div', {}, el('h4', {}, 'COI / BOLD target check'), el('p', {}, supported ? 'Supported at the requested taxon' : human(mito.target_status)), sample.expected_taxon ? el('em', {}, sample.expected_taxon) : null, el('small', {}, 'This checks the supplied identity. It is not an independent deepest-rank taxonomic prediction.')),
    el('div', {}, el('h4', {}, 'Sequence interpretation'), el('p', {}, mito.accepted ? 'Accepted under this run’s policy' : mito.sequence ? 'Candidate retained for review' : 'No sequence recovered'), el('small', {}, 'Circularity, reference-gene recovery and COI agreement are separate evidence. They do not establish every tRNA or exclude every chimera.'))));
  root.append(el('div', {class: 'file-line'}, link(mito.sequence, mito.accepted ? 'Download accepted FASTA ↓' : 'Download candidate FASTA ↓', true), link(mito.annotation.file, 'MITOS annotations'), link(sample.summary_file, 'Evidence summary')));
  if (mito.annotation.features.length) root.append(details(`Annotation table · ${mito.annotation.features.length} feature records`, table(['Gene / feature', 'Type', 'Start', 'End', 'Strand', 'Product'], mito.annotation.features.map(feature => [feature.gene, feature.type, fmt(feature.start), fmt(feature.end), feature.strand, feature.product]))));
  if (mito.hits.length) root.append(details('Supporting COI reference records', table(['BOLD record', 'Reference taxon', 'Identity', 'Aligned bp'], mito.hits.map(hit => [hit.record, hit.taxon, fmt(hit.identity)+'%', fmt(hit.aligned_bp)]))));
  if (mito.files.length) root.append(fileList(mito.files));
  return root;
}
function itsView(sample) {
  const its = sample.its, root = document.createDocumentFragment();
  const index = Math.min(state.locus, Math.max(0, its.loci.length-1)), locus = its.loci[index];
  root.append(metrics([
    ['Reported complete loci', fmt(its.count), '', human(its.status)],
    ['Selected ITS cassette', locus ? fmt(locus.its_length_bp) : '—', 'bp', 'ITS1 · 5.8S · ITS2'],
    ['Full rDNA locus', locus ? fmt(locus.locus_length_bp) : '—', 'bp', 'Includes the reported rRNA flanks'],
    ['Consensus resolution', locus?.taxonomy_rank ? human(locus.taxonomy_rank) : 'Unknown', '', 'From native ITSME consensus evidence']
  ]));
  if (its.loci.length) {
    const select = el('select', {'aria-label': 'Reported ITS locus', onchange: event => {state.locus = Number(event.target.value); renderView();}}, its.loci.map((entry, i) => el('option', {value: i, selected: i === index ? '' : null}, entry.locus)));
    root.append(el('label', {class: 'locus-select'}, 'Reported locus', select));
  }
  root.append(el('div', {class: 'section-head'}, el('div', {}, el('h3', {}, state.mode === 'graph' ? 'The rDNA assembly graph' : 'The reconstructed rDNA locus'), el('p', {}, state.mode === 'graph' ? 'The graph may contain several loci or alternative paths.' : 'Flanking rRNA anchors and the intervening ITS cassette.')), viewSwitch([['map', 'Locus map'], ['graph', 'Assembly graph']])));
  const host = el('div'); root.append(host);
  if (state.mode === 'graph') renderGraph(its.graphs, host);
  else if (locus) host.append(mapSequence(locus.locus_length_bp || Math.max(0, ...locus.regions.map(region => region.end)), locus.regions, true));
  else host.append(empty(its.status === 'not_requested' ? 'ITSME not requested' : ['missing_report','unsupported_report','failed_no_report','not_run'].includes(its.status) ? 'ITSME results unavailable' : 'No complete ITS locus reported', its.status === 'not_requested' ? 'ITSME was not requested for this sample.' : 'Check the run record and retained ITSME files. Partial sequences may still be present.'));
  if (locus) {
    root.append(el('div', {class: 'evidence'}, el('div', {}, el('h4', {}, 'Predicted consensus taxon'), el('em', {}, locus.taxon || 'Unresolved'), el('small', {}, `Resolution: ${locus.taxonomy_rank || 'not recorded'} · ${human(locus.taxonomy_status)}`)), el('div', {}, el('h4', {}, 'Reported cassette lengths'), el('p', {}, `ITS1 ${fmt(locus.its1_length_bp)} bp · 5.8S ${fmt(locus.five_8s_length_bp)} bp · ITS2 ${fmt(locus.its2_length_bp)} bp`), el('small', {}, 'Reported complete loci are kept distinct. Missing assignments are not replaced with the expected sample identity.'))));
    if (locus.consensus_taxonomy) root.append(el('p', {class: 'lineage'}, el('strong', {}, 'Consensus lineage  '), locus.consensus_taxonomy));
  }
  root.append(el('div', {class: 'file-line'}, link(its.sequence, its.sequence?.name.includes('complete_ITS') ? 'ITS cassette FASTA ↓' : 'rDNA loci FASTA ↓', true), its.cassette_sequence?.href !== its.sequence?.href ? link(its.cassette_sequence, 'ITS cassette FASTA') : null, link(its.summary_file, 'ITSME results table')));
  if (its.loci.length) root.append(details('All reported complete loci · '+its.loci.length, table(['Locus', 'ITS bp', 'Full locus bp', 'Consensus taxon', 'Rank'], its.loci.map(entry => [entry.locus, fmt(entry.its_length_bp), fmt(entry.locus_length_bp), entry.taxon, entry.taxonomy_rank]))));
  if (its.files.length) root.append(fileList(its.files)); return root;
}
function runView(sample) {
  const root = document.createDocumentFragment(), qc = sample.qc;
  root.append(metrics([['Input reads', fmt(qc.reads_before), '', 'fastp total reads; not read pairs'], ['Retained reads', fmt(qc.reads_after), '', exists(qc.reads_before) && qc.reads_before > 0 && exists(qc.reads_after) ? fmt(qc.reads_after/qc.reads_before*100)+'% retained' : 'QC report may be unavailable'], ['Q30 after QC', exists(qc.q30_after) ? fmt(qc.q30_after*100)+'%' : '—', '', 'Fraction of retained bases'], ['Recorded stages', fmt(sample.stages.length), '', human(sample.pipeline)+' pipeline']]));
  root.append(el('p', {class: 'run-intro'}, 'A successful tool invocation is separate from biological recovery. Failed or timed-out stages remain visible even when another method produced an accepted sequence.'));
  root.append(el('div', {class: 'progress-track', 'aria-hidden': 'true'}, sample.stages.map(stage => el('div', {class: stage.status.toLowerCase(), title: stage.stage+': '+stage.status}))));
  if (sample.stages.length) root.append(table(['Stage', 'Status', 'Time', 'Note'], sample.stages.map(stage => [human(stage.stage), el('span', {class: 'stage-status '+stage.status.toLowerCase()}, human(stage.status)), numberTime(stage.elapsed_seconds), stage.reason])));
  else root.append(empty('No stage table available', 'Use the saved sample summary and native logs below.'));
  sample.warnings.forEach(message => root.append(el('p', {class: 'note warning'}, message)));
  root.append(el('div', {class: 'file-line'}, link(qc.report, 'Open QC report'), link(sample.summary_file, 'Sample summary')));
  if (sample.logs.length) root.append(fileList(sample.logs, 'Native logs'));
  return root;
}
function numberTime(value) {const n = Number(value); return !exists(value) || !Number.isFinite(n) ? '—' : n < 60 ? fmt(n)+' s' : fmt(n/60)+' min';}
function renderView() {
  const sample = report.samples[state.sample];
  $('view').replaceChildren(state.branch === 'mito' ? mitoView(sample) : state.branch === 'its' ? itsView(sample) : runView(sample));
  if (sample.warnings.length) $('view').prepend(el('p', {class: 'note warning'}, sample.warnings.join(' ')));
  document.querySelectorAll('[data-branch]').forEach(button => {if (button.dataset.branch === state.branch) button.setAttribute('aria-current', 'page'); else button.removeAttribute('aria-current');});
}
function render() {
  const sample = report.samples[state.sample]; sampleList();
  $('sample-name').textContent = sample.id;
  $('sample-position').textContent = `SAMPLE ${String(state.sample+1).padStart(2, '0')} / ${String(report.samples.length).padStart(2, '0')}`;
  $('expected-name').replaceChildren(sample.expected_taxon ? 'Expected identity · ' : 'Expected identity not supplied', sample.expected_taxon ? el('em', {}, sample.expected_taxon) : '');
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
