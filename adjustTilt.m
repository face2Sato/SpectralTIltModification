function [segment, tilt_m] = adjustTilt(segment, tilt_l, rms_p, method, emp_coeff, maxRelErr, fs) 
% スペクトル傾斜を変更するフィルタfractionalSTfilter()はalpha∈[-1,1]を与えられると、非線形にスペクトルを変化させる。
% その為、最初にフィルタが変更可能なスペクトル傾斜の極小値を黄金分割探索で求め、単調増加または単調減少する範囲で最適なalphaを探索する。

% The filter fractionalSTfilter(), which changes the spectral slope, changes the spectrum in a non-linear manner given alpha∈[-1,1].
% Therefore, the first step is to find the minimum value of the spectral slope that can be changed by the filter using a golden split search, 
% and then search for the optimal alpha in the monotonically increasing or monotonically decreasing range.


    % Maximum number of iterations of applying filter
    % If relative error does not become less than maxRelErr no matter how many times the filter is applied, the iteration will be stopped.
    n = 3;

    fc= 500;
    r = 6/5;
    N = floor(log(fs/fc * 0.5) / log(r)); 

    modConfig = struct('fc',fc,'r',r,'N',N);

    alpha_L = -1;
    alpha_R = 1;

    tilt_m = -1000;

    rms_p = rms(segment);

    % While relative error >= maxRelErr, iterate applying filter

    while (abs( (tilt_m - tilt_l) / tilt_l) >= maxRelErr)   &&   n > 0
    
        % Check the limit(tilt_min, tilt_max) of the modification
        [tilt_min, alpha_min] = goldenSectionSearch(segment, alpha_L, alpha_R, rms_p, emp_coeff, method, modConfig, fs);
    
        seg_L = getModifiedSegment(segment, alpha_L, rms_p, emp_coeff, modConfig, fs);
        tilt_L = myGetSpectralTilt(seg_L, fs, method);

        seg_R = getModifiedSegment(segment, alpha_R, rms_p, emp_coeff, modConfig, fs);
        tilt_R = myGetSpectralTilt(seg_R, fs, method);


        if tilt_L < tilt_R
            alpha_max = alpha_R;
            tilt_max = tilt_R;
        else
            alpha_max = alpha_L;
            tilt_max = tilt_L;
        end

        % If lombard tilt exceeds [tilt_min, tilt_max]
        if tilt_l <= tilt_min
 
            segment = getModifiedSegment(segment, alpha_min, rms_p, emp_coeff, modConfig, fs);
            tilt_m = myGetSpectralTilt(segment, fs, method);

        elseif tilt_l >= tilt_max

            segment = getModifiedSegment(segment, alpha_max, rms_p, emp_coeff, modConfig, fs);
            tilt_m = myGetSpectralTilt(segment, fs, method);
        
        else

        % Modify the segment in [-1, alpha_min] or [alpha_min, 1] by Binary search

            if alpha_min < alpha_max
                [tilt_m, segment] = binarySearch(segment, tilt_l, alpha_min, alpha_max, rms_p, emp_coeff, method, modConfig, fs, maxRelErr, 'pos');

            else
                [tilt_m, segment] = binarySearch(segment, tilt_l, alpha_max, alpha_min, rms_p, emp_coeff, method, modConfig, fs, maxRelErr, 'neg');
            end

        end

        n = n-1;
    end

end



function [tilt_extreme, alpha_extreme] = goldenSectionSearch(seg, alpha_l, alpha_r, rms_p, emp_coeff, method, modConfig, fs)

    goldenRatio = (sqrt(5) + 1)/2;

    alpha_m1 = alpha_r - (alpha_r - alpha_l) / goldenRatio;
    alpha_m2 = alpha_l + (alpha_r - alpha_l) / goldenRatio;
    
    seg_m1 = getModifiedSegment(seg, alpha_m1, rms_p, emp_coeff, modConfig, fs);
    tilt_m1 = myGetSpectralTilt(seg_m1, fs, method);

    seg_m2 = getModifiedSegment(seg, alpha_m2, rms_p, emp_coeff, modConfig, fs);
    tilt_m2 = myGetSpectralTilt(seg_m2, fs, method);


    while   alpha_r - alpha_l > eps
        
        if tilt_m1 < tilt_m2
            
            alpha_r = alpha_m2;
            
            alpha_m2 = alpha_m1;
            tilt_m2 = tilt_m1;
            
            alpha_m1 = alpha_r - (alpha_r - alpha_l) / goldenRatio;
        
            seg_m1 = getModifiedSegment(seg, alpha_m1, rms_p, emp_coeff, modConfig, fs);
            tilt_m1 = myGetSpectralTilt(seg_m1, fs, method);
        
        else

            alpha_l = alpha_m1;
            
            alpha_m1 = alpha_m2;
            tilt_m1 = tilt_m2;
            
            alpha_m2 = alpha_l + (alpha_r - alpha_l) / goldenRatio;
        
            seg_m2 = getModifiedSegment(seg, alpha_m2, rms_p, emp_coeff, modConfig, fs);
            tilt_m2 = myGetSpectralTilt(seg_m2, fs, method);
        
        end

    end

    alpha_extreme = alpha_r;

    seg_r = getModifiedSegment(seg, alpha_extreme, rms_p, emp_coeff, modConfig ,fs);
    
    tilt_extreme = myGetSpectralTilt(seg_r, fs, method);

end




function [tilt_mid, seg_mid] = binarySearch(seg, tilt_l, alpha_min, alpha_max, rms_p, emp_coeff, method, modConfig, fs, maxRelErr, pos)
        
    tilt_mid = -1000;

    while abs( (tilt_mid - tilt_l) / tilt_l) >= maxRelErr 
        alpha_mid = (alpha_min + alpha_max) / 2;
        seg_mid = getModifiedSegment(seg, alpha_mid, rms_p, emp_coeff, modConfig, fs);
        tilt_mid = myGetSpectralTilt(seg_mid, fs, method);
    
        if strcmp(pos,'pos')
            
            if tilt_mid < tilt_l
                alpha_min = alpha_mid;
            else
                alpha_max = alpha_mid;
            end

        else

            if tilt_mid < tilt_l
                alpha_max = alpha_mid;
            else
                alpha_min = alpha_mid;
            end

        end
    end

end



function seg_mod = getModifiedSegment(seg, alpha, rms_p, emp_coeff, modConfig, fs) 

    % Remove DC offset
    seg = seg - mean(seg);

    % pre-emp
    seg = filter([1 emp_coeff], 1, seg);

    seg_mod = fractionalSTfilter(seg, fs, alpha, modConfig.fc, modConfig.r, modConfig.N);

    %de-emp
    seg_mod = filter(1, [1 emp_coeff], seg_mod);

    % gain match to lombard segment
    seg_mod = seg_mod * rms_p/rms(seg_mod);

end