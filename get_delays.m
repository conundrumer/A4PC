function [ delays_out, iter, blur ] = get_delays( A, B, AL, AR, BL, BR, d_max )

DEBUG = false;

RATIO = 8;
MAX_JUMP = 3*RATIO/2;
MAX_CURVE = 3*RATIO/2;
MIN_DISCONTINUITY_JUMP = 64*RATIO/2; % 64 samples
N = 256;
HOP = 256;
RANGE = 4;
BLUR_RANGE = 16;
smooth_kernel = fspecial('gaussian', 7, 0.5);
h0 = N*RATIO*2^(RANGE-1)/4;
num_frames = floor(length(A) / HOP) + 1;
% d_min = 3*RATIO/2;

big_jump = false;

function [M_, ds_] = get_peaks(c_filtered)
  [M_ ds_] = max(c_filtered,[],1);
  M_ = M_(2:end-1);
  ds_ = ds_(2:end-1);
end

function [jumps, spikes] = get_discontinuities(delays)
  jumps = abs(diff(delays)) > MAX_JUMP;
  spikes = abs(diff(diff(delays))) > MAX_CURVE;
end

function [is_smooth] = smooth_enough(jumps, spikes, delays)
  num_jumps = sum(jumps);
  num_spikes = sum(spikes);
  js = abs(diff(delays));

  if (num_jumps == 1 & num_spikes == 2)
    jump_loc = find(jumps);
    spike_loc = find(spikes, 2);
    if (spike_loc(1) + 1 == spike_loc(2) & jump_loc == spike_loc(2))
      % this is a clean discontinuity
      js(jump_loc);
      is_smooth = (js(jump_loc) >= MIN_DISCONTINUITY_JUMP);
      if (is_smooth)
        big_jump = jump_loc;
        return;
      end
    end
  end

  is_smooth = ((num_spikes < 2 & num_jumps < 1)); % no spikes/jumps or 1 big jump
  % is_smooth = ((num_spikes < 2 & num_jumps < 1) | (sum(js(jumps)) < MIN_DISCONTINUITY_JUMP/2)); % no spikes/jumps or 1 big jump
end

c_multi = zeros(N*RATIO/4, floor(length(A) / HOP) + 1);
weights = zeros(h0, 1);

min_breaks = 2*num_frames;

iter = 0;
for i = 1:RANGE
  blur = 0;

  n = N*2^(i-1);
  h = n * RATIO / 4;

  A_padded = [AL(length(AL)-n/2+1:length(AL)) ; A ; AR(1:n/2)];
  B_padded = [BL(length(BL)-n/2+1:length(BL)) ; B ; BR(1:n/2)];

  % A_padded = padarray(A, n/2);
  % B_padded = padarray(B, n/2);

  c = get_interference_cepstrum(A_padded, B_padded, n*1, HOP, RATIO);

  prev_h = size(c_multi, 1);

  c(1:prev_h,:) = c_multi + c(1:prev_h,:);

  c_multi = c;

  weights(1:h) = weights(1:h) + ones(h,1);

  c_normed = bsxfun(@rdivide, c_multi, weights(1:h));
  % c_filtered = c_normed(d_min+1:min(h,RATIO/2*d_max),:); % limit possible delays to reasonable ones
  c_filtered = c_normed(1:min(h,RATIO/2*d_max),:); % limit possible delays to reasonable ones
  [M, delays] = get_peaks(c_filtered);
  % delays = delays + d_min;

  % subplot(RANGE, 2, 2*(i-1)+1);
  % plot(delays);
  % subplot(RANGE, 2, 2*(i-1)+2);

  [jumps, spikes] = get_discontinuities(delays);

  if (smooth_enough(jumps, spikes, delays)) % only allow one discontinuity
    iter = i;
    delays_out = delays;
    c_best = c_filtered;
    break;
  end

  if (sum(jumps) > num_frames/4 & (i < RANGE)) % a lot => can't save through blurring
    continue;
  end

  done = false;
  for j = 1:BLUR_RANGE % blur until
    blur = j;
    c_filtered = imfilter(c_filtered, smooth_kernel);
    [M, delays] = get_peaks(c_filtered);

    [jumps, spikes] = get_discontinuities(delays);

    if (smooth_enough(jumps, spikes, delays)) % only allow one discontinuity
      done = true;
      break;
    end
  end

  if (done)
    iter = i;
    delays_out = delays;
    c_best = c_filtered;
    break;
  end

  num_breaks = sum(jumps) + sum(spikes);

  if ((num_breaks < min_breaks) | (i == RANGE));
    iter = i;
    min_breaks = num_breaks;
    c_best = c_filtered;
    delays_out = delays;
  end
end

% stem(jumps.*diff(delays));
% hold on;
% % stem(spikes.*diff(diff(delays)));
% stem([1.5:1:length(spikes)+0.5],spikes.*diff(diff(delays)));
% hold off;

% c_normed = bsxfun(@rdivide, c_multi, weights);

% c_smoothed = imfilter(c_normed, smooth_kernel);

% [M, delays] = max(c_smoothed,[],1);

delays_out = delays_out';
M = M';

function [ds_out, rs, ms] = smooth_delays(ds_in, ms)

  ds = ds_in;
  % ms = ms';

  avg = sum(ms.*ds)./sum(ms);
  dev = sqrt(var(ds,ms));


  keep = ( abs(ds - avg) < 2.5*dev );

  ds = round(naninterp(ds .* keep)); % remove by turning to zero

  ms = ms .* keep; % remove all weight

  range = [1:length(ds)]';
  coeffs = lscov([ones(size(range)) (range-1) (range-1).^2 ], ds, ms.^2);
  fit = coeffs(1) + coeffs(2)*(range-1) + coeffs(3)*(range-1).^2;

  ds = ds + fit .* (ds == 0);

  ds_out = fit;

  rs = ds - fit;

end

function [ds_out, rs_out] = add_flutter(ds, rs, ms)

  ds_out = ds;
  rs_out = rs;

  k = 2;
  window = gausswin(length(rs));
  pxx = 10*log10(periodogram(rs, window, length(rs) * k));
  [peaks peaks_i peaks_w peaks_p] = findpeaks(pxx(15*k: 25*k));
  [peak max_peaks_i] = max(peaks);
  peak_i = peaks_i(max_peaks_i);
  peak_p = peaks_p(max_peaks_i);

  peak_i = peak_i + 15*k - 1;
  avg = mean(pxx(1: 25*k));

  flutter = [];

  thresh = 10 + 10*log10(6 / RATIO);

  if DEBUG
    figure(1);
    subplot(2,1,2);
    plot(pxx(1:25*k));

    hold on;
    stem(peak_i, peak);
    line([1 25*k], [avg avg]);
    hold off;
  end

  if (peak > thresh & peak_p > thresh & peak - avg > (thresh - 3))
    y = pxx(peak_i-1:peak_i+1);
    xincrement = 0.5 * (y(1) - y(3)) / (y(1) - 2 * y(2) + y(3));

    flutter_freq = (xincrement + peak_i)/2 - 1;

    % t = [0:length(ds)-1]';
    % waves = [cos(flutter_freq * t) sin(flutter_freq * t)];

    % % coeffs = lscov(waves, ds, M .* window);
    % coeffs = lscov(waves, ds, M);
    % flutter = (waves * coeffs)';
    filter_len = 32;
    bandwidth = 4;
    % res = 0.75;

    freq_range_band = [ flutter_freq - bandwidth/2, flutter_freq + bandwidth/2 ];
    H_band = fir1(filter_len, freq_range_band / (length(rs)/2));

    freq_range_low = flutter_freq/3;
    H_low = fir1(filter_len, freq_range_low / (length(rs)/2));

    % H = (H_band*res + H_low*(1-res));
    H = H_band + H_low;

    flutter = filter(H, 1, [rs; zeros(filter_len, 1)]);
    flutter = flutter(filter_len/2+1:end-filter_len/2);

    % weighted average between flutter and residue
    w = ms / max(ms);
    w = w.^2;
    flutter = (flutter.*(1-w)) + (rs.*w);

    % filter again to make smoother
    flutter = filter(H, 1, [rs; zeros(filter_len, 1)]);
    flutter = flutter(filter_len/2+1:end-filter_len/2);

    if DEBUG
      hold on;
      stem(xincrement + peak_i, peak);
      hold off;
      % figure(2);
      % plot(residues);
      % hold on;
      % plot(flutter);
      % hold off;
    end

    ds_out = ds+flutter;
    rs_out = rs - flutter;
  end

end

if big_jump
  [left_delays left_residues left_m] = smooth_delays(delays_out(1:big_jump), M(1:big_jump));
  [right_delays right_residues right_m] = smooth_delays(delays_out(big_jump+1:end), M(big_jump+1:end));
  delays_out = [left_delays; right_delays];
  residues = [left_residues; right_residues];
  M = [left_m; right_m];
else
  [delays_out residues M] = smooth_delays(delays_out, M);
end

[delays_out residues] = add_flutter(delays_out, residues, M);


if DEBUG
  figure(1);

  delays_old = delays_out;

  subplot(2,1,1);

% if 0
  c_norm = c_best;
  % c_norm = bsxfun(@rdivide, c_best, M);

  avg = round(mean(delays_out));
  span = 20;
  % span = max(10,round(2*std(delays_out)));

  low = max(1, round(min(delays_out))-span);
  high = min(size(c_norm, 1), round(max(delays_out))+span);

  axis ij;
  imagesc(c_norm(low:high,2:end-1));
  title(['iter: ' int2str(iter) ', blur: ' int2str(blur)]);

  % range = [1:num_frames]';
  % coeffs = lscov([ones(size(range)) (range-1) (range-1).^2 ], delays_out', M);
  % fit = coeffs(1) + coeffs(2)*(range-1) + coeffs(3)*(range-1).^2;

  % delays_out = delays_out + fit' .* (delays_out == 0);

  hold on;
  % plot(delays_old-low+1, 'color', 'green');
  plot(delays_out+residues -low+1, 'color', 'red');

  plot(delays_out -low+1, '--', 'color', 'black');

  % if (length(flutter) > 0)
  %   plot(delays_out+flutter -low+1, 'color', 'black');

  %   delays_out = delays_out+flutter;
  %   residues = residues - flutter;
  % end
  % plot(fitter-low+1, '--', 'color', 'black');
  % plot(fitish-low+1, '--', 'color', 'black');
  hold off;
  set(gca,'YDir','normal')

  % title([ 'coefficients: ', num2str(coeffs(1)), ' ', num2str(coeffs(2)), ' ', num2str(coeffs(3)) ]);

  drawnow;
  waitforbuttonpress;
end

delays_out = max(0, (delays_out-1) / (RATIO/2));

end

