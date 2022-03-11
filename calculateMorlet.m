function morletWave = calculateMorlet(time,f)
sigma = 4.5/(2*pi*f);
%morletWave = exp(2*pi*1i*10.*time) .* exp(-time.^2./(2*sigma^2));
morletWave = cmorwavf(round(-length(time)/2),round(length(time)/2),length(time),f,2*sigma);