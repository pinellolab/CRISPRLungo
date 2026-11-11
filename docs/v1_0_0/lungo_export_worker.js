/* Runs minimap2 for lungo_export.js (CLI-identical insertion annotation of inserts >= 100 bp). */
importScripts('https://biowasm.com/cdn/v3/aioli.js');

self.onmessage = async function (event) {
  try {
    var CLI = await new Aioli(['minimap2/2.22'], { printInterleaved: false });
    await CLI.mount([{ name: 'reference_for_ins.fasta', data: event.data.reference }]);
    await CLI.mount([{ name: 'large_insertion.fasta', data: event.data.query }]);
    var out = await CLI.exec('minimap2 -ax map-ont reference_for_ins.fasta large_insertion.fasta');
    self.postMessage({ ok: true, sam: out.stdout });
  } catch (err) {
    self.postMessage({ ok: false, error: String(err && err.message || err) });
  }
};
