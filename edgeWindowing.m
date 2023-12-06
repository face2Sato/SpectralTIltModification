function output = edgeWindowing(input, windowsizeSample, option)
% Do windowing the edges of the input signal
% inputdata        : input signal
% windowsizeSample : Specify window size by number of samples (must be even number)

    if nargin < 3
        option = 'both';
    end

    window = hann(windowsizeSample,'periodic');

    if strcmp(option,'left')

        window_bef = window(1:(windowsizeSample/2));
        newwindow = [window_bef ; ones((length(input)-windowsizeSample/2),1)];
        
    elseif strcmp(option,'right')

        window_aft = window((windowsizeSample/2+1):windowsizeSample);
        newwindow = [ones((length(input)-windowsizeSample/2),1) ; window_aft];
        
    elseif strcmp(option,'inverse')
    
        window_bef = window((windowsizeSample/2+1):windowsizeSample);
        window_aft = window(1:(windowsizeSample/2));
        newwindow = [window_bef ; zeros((length(input)-windowsizeSample),1) ; window_aft];
    
    else
    
        window_bef = window(1:(windowsizeSample/2));
        window_aft = window((windowsizeSample/2+1):windowsizeSample);
        newwindow = [window_bef ; ones((length(input)-windowsizeSample),1) ; window_aft];
        

    end 

    output = input .* newwindow;

end