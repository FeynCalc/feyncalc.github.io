#!/usr/bin/env node
// Search indexes for the manuals and the example gallery, and the index of the
// ⌘K palette.
//
//     node tools/search-index.js        (from the repository root)
//
// Run it after the manuals or the examples were regenerated and commit what it
// writes. Needs Node 16 or newer and nothing else: the engine and the HTML
// splitter come from the sitesearch submodule
// (git submodule update --init).
//
// Why a script and not Jekyll: the manuals are static HTML made by pandoc.
// Jekyll copies them as they are and Liquid cannot read their text, so their
// index is made here and lies ready as a file. The site's own Jekyll pages are
// indexed by Jekyll itself (sitesearch/jekyll/search-index.json, which serves
// /search-index.json) and that index names these files as its `shards` —
// see _data/sitesearch.yml.
//
// What it writes:
//   <book>/search-index.json, <book>/search-index-text/   the text, in two tiers
//   nav-index.json                                         page names, for ⌘K

const fs = require('fs');
const path = require('path');
const { body, blocks, text } = require('../sitesearch/node/html.js');
const { split } = require('../sitesearch/node/two-tier.js');

const ROOT = path.resolve(__dirname, '..');

// The parts of the site that are searched, in the order the results appear.
// `home` is the site's page for the section: the palette lands there when the
// section itself is asked for.
const BOOKS = [
	{ dir: 'FeynCalcBook', section: 'FeynCalc manual', order: 1, home: '/reference', homeTitle: 'FC Guide (stable)' },
	{ dir: 'FeynCalcBookDev', section: 'FeynCalc manual (development)', order: 2, home: '/referenceDev', homeTitle: 'FC Guide (development)' },
	{ dir: 'FeynHelpersBookDev', section: 'FeynHelpers manual (development)', order: 3, home: '/referenceFeynHelpersDev', homeTitle: 'FH Guide (development)' },
	{ dir: 'FeynCalcExamples', section: 'Examples', order: 4, home: '/examples', homeTitle: 'Example gallery' },
];

// Site pages for the palette. They are few and named by hand: the rest of the
// site is the manuals.
const SITE = [
	{ url: '/', title: 'Home' },
	{ url: '/search', title: 'Search' },
	{ url: 'https://github.com/FeynCalc/feyncalc/discussions', title: 'Forum' },
	{ url: '/OldForum/thread.html', title: 'Old forum archive' },
	{ url: '/impressum', title: 'Imprint (Impressum)' },
	{ url: '/datenschutz', title: 'Privacy (Datenschutz)' },
];

// The text of the search index is shown and searched as it stands, so the
// manual's generic heading ("FeynCalc manual (development version)" on every
// page) goes. Formulas stay: html.js turns pandoc's math spans into \(…\) and
// \[…\], the engine matches them however they are spelled (γ^μγ^ν finds
// \bar{\gamma}^{\mu } \bar{\gamma}^{\nu }) and draws them with KaTeX. Code
// stays too — the examples are how a function is used.
const FURNITURE = [
	/<header id="title-block-header">[\s\S]*?<\/header>/g,
	/<p>\$\$\$\$<\/p>/g,
];

// Except the largest output cells. Past 600 characters the engine shows a
// formula as cut TeX rather than drawing it, and nobody types a query that only
// a three-screen expression answers; the 232 displayed formulas over 1000
// characters (of 12 000) are a quarter of all the TeX in the manuals.
const HUGE_MATH = 1000;
const DISPLAY = /<span\s+class="math display">([\s\S]*?)<\/span>/g;
const dropHuge = (html) => html.replace(DISPLAY, (m, tex) => (tex.length > HUGE_MATH ? '' : m));

// A chunk of text above this size is split into a map of words and the
// chunks behind it (sitesearch/FORMAT.md, "two tiers").
const CHUNK = 32;

function walk(dir) {
	const out = [];
	for (const e of fs.readdirSync(dir, { withFileTypes: true })) {
		const p = path.join(dir, e.name);
		if (e.isDirectory()) {
			if (!['css', 'js', 'img', 'search-index-text'].includes(e.name)) out.push(...walk(p));
		} else if (e.name.endsWith('.html')) {
			out.push(p);
		}
	}
	return out;
}

const byName = new Intl.Collator('en', { numeric: true }).compare;

// A page's name: the example's own description where it states one, else the
// first heading, else the file name — which in a manual is the symbol itself.
function nameOf(html, rel) {
	const d = html.match(/description <span class="ex">=<\/span> <span class="st">&quot;([\s\S]*?)&quot;<\/span>/);
	if (d) return text(d[1]);
	const h = html.match(/<h2\b[^>]*>([\s\S]*?)<\/h2>/);
	if (h) return text(h[1]);
	return path.basename(rel, '.html');
}

function build(book) {
	const dir = path.join(ROOT, book.dir);
	const files = walk(dir)
		.map((f) => path.relative(ROOT, f).split(path.sep).join('/'))
		// Overview and guides first, then the symbols in alphabetical order.
		.sort((a, b) => (b.includes('/Extra/') - a.includes('/Extra/')) || byName(a, b));

	const pages = files.map((rel) => {
		const html = fs.readFileSync(path.join(ROOT, rel), 'utf8');
		const url = '/' + rel;
		const title = nameOf(html, rel);
		const base = path.basename(rel, '.html');
		// Text results come after every name that matched (the search page
		// reads nav-index.json first, at the orders below ten).
		const page = { url, title, order: 10 + book.order, blocks: [] };
		if (base !== title) page.also = base;
		for (const b of blocks(body(dropHuge(html), FURNITURE))) {
			if (b.text.length < 3) continue;
			const block = { text: b.text };
			if (b.anchor) block.anchor = b.anchor;
			if (b.section && b.section !== title) block.section = b.section;
			page.blocks.push(Object.keys(block).length === 1 ? block.text : block);
		}
		return page;
	});

	const index = { pages };
	const flat = JSON.stringify(index).length;
	const out = split(index, { chunk: CHUNK, repeats: 4, text: 'search-index-text/' });

	const textDir = path.join(dir, 'search-index-text');
	fs.rmSync(textDir, { recursive: true, force: true });
	fs.mkdirSync(textDir);
	let textBytes = 0;
	out.chunks.forEach((chunk, id) => {
		const s = JSON.stringify(chunk);
		textBytes += s.length;
		fs.writeFileSync(path.join(textDir, id + '.json'), s);
	});
	const map = JSON.stringify(out.map);
	fs.writeFileSync(path.join(dir, 'search-index.json'), map);

	const kb = (n) => Math.round(n / 1024) + ' KB';
	console.log(`/${book.dir}/search-index.json: ${pages.length} pages, ${out.words} words; ` +
		`flat ${kb(flat)} -> map ${kb(map.length)} + text ${kb(textBytes)} in ${out.chunks.length} chunks`);

	return pages.map((p) => ({ url: p.url, title: p.title, section: book.section, order: book.order }));
}

const nav = SITE.map((p) => ({ ...p, section: 'FeynCalc', order: 0 }));
for (const book of BOOKS) {
	// The site's page for the section is its front door; the manual's own
	// pages sit one level deeper and all at the same depth, so without it every
	// one of them would be taken for the section's landing page.
	nav.push({ url: book.home, title: book.homeTitle, section: book.section, order: book.order });
	for (const p of build(book)) nav.push(p);
}
const navJson = JSON.stringify({ pages: nav });
fs.writeFileSync(path.join(ROOT, 'nav-index.json'), navJson);
console.log(`/nav-index.json: ${nav.length} pages, ${Math.round(navJson.length / 1024)} KB`);
