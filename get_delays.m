function [ delays ] = get_delays( A, B )

N = 256;
HOP = 256;
RATIO = 6;
RANGE = 6;
smooth_kernel = fspecial('gaussian', [16 16], 3);
h0 = N*RATIO*2^(RANGE-1)/4;

c_multi = zeros(h0, floor(length(A) / HOP) + 1);
weights = zeros(h0, 1);

for i = 1:RANGE

  n = N*2^(i-1);
  h = n * RATIO / 4;

  c = get_interference_cepstrum(padarray(A, n/2), padarray(B, n/2), n*1, HOP, RATIO);

  c_multi(1:h,:) = c_multi(1:h,:) + c;

  weights(1:h) = weights(1:h) + ones(h,1);
end

c_normed = bsxfun(@rdivide, c_multi, weights);

c_smoothed = imfilter(c_normed, smooth_kernel);

[ms, delays] = max(c_smoothed,[],1);

% debug
% if 1
if 0
  c_norm = bsxfun(@rdivide, c_smoothed, ms);

  avg = round(mean(delays));
  range = 10;
  % range = max(10,round(2*std(delays)));

  % low = max(1, avg-range);
  % high = min(h0, avg+range);

  delays_trunc = delays(2:length(delays)-1);

  low = max(1, min(delays_trunc)-range);
  high = min(h0, max(delays_trunc)+range);

  axis ij;
  imagesc(c_norm(low:high,:));

  hold on;
  plot(delays-low+1, 'color', 'red');
  hold off;
  set(gca,'YDir','normal')
end

delays / 3;

end

