import init, { analyze } from './pkg_v5/CRISPRlungo_regular.js?v=1.1.0';

// 20261111_cachebust: the wasm is cached hard by browsers, so a rebuilt engine can
// keep running the old binary until someone clears the cache. Tying the URL to the
// app version means a version bump fetches the new one. Keep this in step with
// CRISPRlungo_log.VERSION and the navbar badge.
const ENGINE_VERSION = '1.1.0';

self.onmessage = async (e) => {

    let {
        controlSamString,
        treatedSamString,
        inducedMutationString,
        reference,
        cvPos,
        cvPos2,
        window,
        wholeWindow,
        filter1,
        windowFilter,
        alpha,          // 20261111_worker: these were sent but dropped on the floor
        length_min,
    } = e.data;

    if (cvPos2 == "") {
        cvPos2 = undefined;
    }

    try {

        await init('./pkg_v5/CRISPRlungo_regular_bg.wasm?v=' + ENGINE_VERSION);

        const result = analyze(
            controlSamString,
            treatedSamString,
            inducedMutationString,
            reference,
            cvPos,
            cvPos2,
            window,
            wholeWindow,
            filter1,
            windowFilter,
            (typeof alpha === 'number' ? alpha : 0.05),      // 20261111_worker
            (typeof length_min === 'number' ? length_min : 10)  // 20261111_worker
        );


        self.postMessage({ type: 'success', result });
    } catch (error) {
        // 20261111_worker: 'error_run' was the old name and nothing checked for it,
        // so a wasm panic came back looking like a successful run with no result.
        self.postMessage({ type: 'error', kind: 'error_run',
                           message: error.message, stack: error.stack });
    }
};
