% input:
% output: v_uv...0(voiceless),1(voiced)

% wakaran

function  = detectVUV(frame, fs, threshold)
    method = 6;
    emp_coeff = -0.97;
    weight = 0.8;

    spectilt = myGetSpectralTilt(frame, fs, method, emp_coeff);
    intensity = rms(frame);

    voiced = intensity - weight * spectilt;
end