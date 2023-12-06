function seg_m = LPCTransplantation(seg_p, seg_l, emp_coeff, fs)
%     lpcOrder = fix(fs/1000) + 5;
    lpcOrder = 13;

    rms_p = rms(seg_p);

    seg_p = seg_p - mean(seg_p);
    seg_l = seg_l - mean(seg_l);

    seg_p = filter([1 emp_coeff], 1, seg_p);
    seg_l = filter([1 emp_coeff], 1, seg_l);

    ap = lpc(seg_p, lpcOrder);
    al = lpc(seg_l, lpcOrder);

    seg_m = filter([1 ap(2:end)], [1 al(2:end)], seg_p);
%     seg_m = filter(ap, al, seg_p);

    seg_m = filter(1, [1 emp_coeff], seg_m);
    seg_m = seg_m * rms_p / rms(seg_m);

    
end