import init, { analyze } from './pkg_v5/CRISPRlungo_regular.js';

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
    } = e.data;

    if (cvPos2 == "") {
        cvPos2 = undefined;
    }

    try {

        await init('./pkg_v5/CRISPRlungo_regular_bg.wasm');

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
            0.05,
            10
        );


        self.postMessage({ type: 'success', result });
    } catch (error) {
        self.postMessage({ type: 'error_run', message: error.message, stack: error.stack });
    }
};
