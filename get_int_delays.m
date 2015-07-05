function [ delays ] = get_int_delays( A, B, n, hop )
% get_int_delays: given signals A and B ,
% returns the instantaeous delays as integers
% adjusted for sign and with outliers removed and interpolated
RATIO = 2;
MAX_DELAY_DIFF = 2;

cepstrum = get_interference_cepstrum(A, B, n, hop, RATIO);

[M, I] = max(cepstrum,[],1);

delays = 2 * (I - 1) / RATIO;
avg = mean(delays);

max_amp_diffs = []; % for debugging
amp_LR_diffs = []; % for debugging

old_delays = delays;
delays(1) = 0;

m = hop * 1;
% attach sign
prev_i = 1;
avg_d = 0;
for i = 2:length(delays)

  x = hop * (i - 1) + 1; % offset
  d = delays(i);

  if x-d < 1 | x+m+d > length(B)
    delays(i) = NaN; % out of range delay
    continue;
  end

  frame_A = A(x:x+m-1);

  % relative to A...
  frame_BL = B(x-d:x+m-d-1); % negative delay
  frame_BR = B(x+d:x+m+d-1); % positive delay
  if x-d < 1
    frame_BL = zeros(m, 1);
  end
  if x+m+d > length(B)
    frame_BR = zeros(m, 1);
  end

  amp_A = mean(abs(frame_A));
  amp_ABL = mean(abs(frame_A - frame_BL));
  amp_ABR = mean(abs(frame_A - frame_BR));

  % statistically likely that mixes increase in volume by sqrt(2)
  % UNLESS they're in phase by within one sample
  if amp_ABL >= amp_A & amp_ABR >= amp_A & (abs(d - delays(prev_i)) > MAX_DELAY_DIFF & abs(-d - delays(prev_i)) > MAX_DELAY_DIFF)
    delays(i) = NaN; % this delay doesn't do any phase cancellation and is an outlier
    continue;
  end

  if amp_ABL < amp_ABR % negative delay has more effect
    delays(i) = -d;
    max_amp_diffs(i) = amp_A / (amp_A - amp_ABL);
  else
    delays(i) = d;
    max_amp_diffs(i) = amp_A / (amp_A - amp_ABR);
  end

  if amp_ABL ~= amp_ABR
    amp_LR_diffs(i) = amp_ABL / amp_ABR;
  end

  % flip signs if more optimal result
  d = delays(i);
  if abs(-d - delays(prev_i)) < abs(d - delays(prev_i)) & abs(-d - avg_d) <= MAX_DELAY_DIFF
    delays(i) = -d;
  end

  prev_i = i;
  avg_d = (avg_d + delays(i)) / 2;

end

% remove outliers, reinsert ok points
prev_i = 1;
i = 2;
for next_i = 3:length(delays)
  if isnan(delays(next_i))
    continue;
  end

  prev_diff = abs(delays(prev_i) - delays(i));
  next_diff = abs(delays(i) - delays(next_i));

  if prev_diff > MAX_DELAY_DIFF & next_diff > MAX_DELAY_DIFF
    delays(i) = NaN;
  else
    prev_i = i;
  end

  i = next_i;

end

% interpolate missing delays or replace with 0
delays(length(delays)) = 0;
delays = fixgaps(delays);
delays = round(delays);


% plot(old_delays);
% hold on;

% % stem(avg * log(max_amp_diffs), 'x');
% % stem(avg * abs(log(amp_LR_diffs)), 'x');
% plot(delays, '-o');

% hold off;

end

