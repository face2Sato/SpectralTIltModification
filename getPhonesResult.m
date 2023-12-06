function result_phones = getPhonesResult(result)

    phones = unique(replace(result.phone,{'_I','_S','_B','_E'},''));
    
    means = zeros(length(phones),1);
    stds = zeros(length(phones),1);
    
    for i = 1:height(phones)
        means(i) = mean(abs(result{contains(result.phone, phones(i)),'relerr_seg'}));
        stds(i) = std(result{contains(result.phone, phones(i)),'relerr_seg'});
    end
    
    result_phones = table(phones,means,stds);

end