/* DOM-only interaction checks. Requires jsdom for tests, not for report users. */
const fs = require('node:fs');
const assert = require('node:assert/strict');
const {JSDOM, VirtualConsole} = require(process.env.TORTUGA_JSDOM || 'jsdom');
const errors = [];
const consoleCapture = new VirtualConsole();
consoleCapture.on('jsdomError', error => errors.push(error.message));
const dom = new JSDOM(fs.readFileSync(process.argv[2], 'utf8'), {
  runScripts: 'dangerously', url: 'https://offline-report.invalid/', virtualConsole: consoleCapture,
  beforeParse(window) {
    window.HTMLDialogElement.prototype.showModal = function() {this.open = true;};
    window.HTMLDialogElement.prototype.close = function() {this.open = false;};
    window.fetch = () => {throw new Error('An offline report must not fetch resources');};
    window.XMLHttpRequest = class {constructor() {throw new Error('An offline report must not request resources');}};
  }
});
const document = dom.window.document;
const click = selector => {const node = document.querySelector(selector); assert(node, 'Missing control '+selector); node.dispatchEvent(new dom.window.MouseEvent('click', {bubbles:true}));};
const text = () => document.querySelector('#view').textContent;
assert.equal(document.querySelector('#sample-name').textContent, 'W30877');
assert.equal(document.querySelectorAll('.sample-row').length, 8);
assert(text().includes('15,740'));
assert(document.querySelectorAll('.gene').length > 10);
click('.gene');
assert(document.querySelector('.feature-readout').textContent.includes('cox1'));
click('.view-switch button:last-child');
assert.equal(document.querySelectorAll('.graph .segment').length, 12);
click('.graph .segment');
assert(document.querySelector('.graph-readout').textContent.includes('Segment 1'));
click('[aria-label="Zoom in"]'); click('[aria-label="Zoom out"]'); click('[aria-label="Fit graph"]');
click('[data-branch="its"]');
assert(text().includes('990'));
assert.equal(document.querySelectorAll('.gene').length, 5);
const select = document.querySelector('[aria-label="Reported ITS locus"]');
select.value = '1'; select.dispatchEvent(new dom.window.Event('change'));
assert(text().includes('1,040'));
click('.view-switch button:last-child');
assert.equal(document.querySelectorAll('.graph .segment').length, 8);
click('[data-branch="run"]');
assert(text().includes('timeout'));
click('#compare'); assert(document.querySelector('#comparison').open);
click('#close-comparison'); assert(!document.querySelector('#comparison').open);
const input = document.querySelector('#search');
input.value = 'Geodia'; input.dispatchEvent(new dom.window.Event('input'));
assert.equal(document.querySelectorAll('.sample-row').length, 1);
click('.sample-row'); click('[data-branch="its"]');
assert(text().includes('No complete ITS locus reported'));
assert(document.querySelector('#sample-status').textContent.includes('ITSME · completed'));
assert(!document.querySelector('#sample-status').classList.contains('fail'));
click('[data-branch="mito"]');
assert(text().includes('No annotation map available'));
assert(!text().includes('null'));
assert.equal(document.querySelectorAll('[src^="http"],link[href^="http"]').length, 0);
assert.deepEqual(errors, []);
console.log('Passed offline DOM checks: filtering, sample/branch/locus switching, gene selection, GFA/FASTG selection, zoom buttons, comparison, run failures and missing results.');
dom.window.close();
