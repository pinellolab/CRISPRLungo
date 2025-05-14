import init, { analyze } from './pkg/CRISPRlungo_regular.js';

self.onmessage = async (e) => {

    const {
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

    try {

        await init('./pkg/CRISPRlungo_regular_bg.wasm');

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
            windowFilter
        );

        self.postMessage({ type: 'success',  result    });
    } catch (error) {
        self.postMessage({ type: 'error_run', message: error.message, stack: error.stack });
    }
};
