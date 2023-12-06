function seg_m = iterateSTfilter(seg_p, fs, tilt_p, tilt_l, method)

    fc= 450;
    r = 6/5;
    N = ceil(log(fs/fc * 0.5) / log(r));
    [b,a] = butter(1,fc/(fs/2));
    alpha = 1000;
    k = 3;
    rms_p = rms(seg_p);
    seg_p = seg_p - mean(seg_p);
    seg_m = seg_p;
    tilt_m = tilt_p;


    while abs(alpha) >= 1 & k > 0

        alpha = (tilt_l - tilt_m)/20 * log2(10);
        
        if alpha > 1
            seg_m = filter(a,b,seg_m);
        elseif alpha < -1
            seg_m = filter(b,a,seg_m);
        end
        
        seg_m = seg_m - mean(seg_m);
        seg_m = seg_m * rms_p/rms(seg_m);
        
        tilt_m = myGetSpectralTilt(seg_m, fs, method);                   
        k = k-1;
    end

    if alpha > 1, alpha = 1; end
    if alpha < -1, alpha = -1; end
    
    seg_m = fractionalSTfilter(seg_m, fs, alpha, fc, r, N);  
    seg_m = seg_m - mean(seg_m);
    seg_m = seg_m * rms_p/rms(seg_m);
end