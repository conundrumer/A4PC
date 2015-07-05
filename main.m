tic
display('loading files');
tic
[A, fs] = audioread('./audio/full_A.wav');
A = A(:,1) + A(:,2);

[B, fs] = audioread('./audio/full_B.wav');
B = B(:,1) + B(:,2);

A = A(1:length(B));

BIG_N = 2048*8;
BIG_HOP = 2*BIG_N;
HOP = 256;
toc
%%%%%%%%% coarse analysis
display('coarse analyiss');
tic
m = 2 * BIG_HOP / HOP; % number of frames per chunk
[A_ B_ ds] = normalize_delay(A, B, BIG_N, BIG_HOP);
A_ = padarray(A_, 2*BIG_HOP);
B_ = padarray(B_, 2*BIG_HOP);

delay_offset = min(ds);
toc
%%%%%% medium analysis
display('analyzing delays');
tic
window = hann(m, 'periodic');
summed_delays = zeros(m*(length(ds)+4)/2,1);



% for i = 104:106 % turbulence
% for i = 183:184 % jump
% for i = 168:168 % bad jump
% for i = 243:length(ds)-2 % more wiggles
% for i = 175:190
for i = 1:length(ds)
% for i = 1:10
  j = 2*BIG_HOP + (i-1)*BIG_HOP;

  delay_set = [ds(max(1,i-1)), ds(i), ds(min(length(ds), i+1))];

  d = min(delay_set) - round(delay_offset/2);

  d_max = max(delay_set);

  left = j+1-BIG_HOP-HOP;
  right = j+BIG_HOP;

  frame_A = A_(left:right);
  frame_AL = A_(max(1,left-BIG_HOP):j);
  frame_AR = A_(right+1:right+BIG_HOP);

  frame_B = B_(d+left:d+right);
  frame_BL = B_(max(1,d+left-BIG_HOP):d+j);
  frame_BR = B_(d+right+1:d+right+BIG_HOP);

  [delays, iter, blur] = get_delays(frame_A, frame_B, frame_AL, frame_AR, frame_BL, frame_BR, d_max);

  offset = m + (i-1)*m/2;
  L = offset+1-m/2;
  R = offset+m/2;

  summed_delays(L:R) = summed_delays(L:R) + window.*(delays(1:end) + d);

end


toc
% %%%%%%%%  fine analysis
display('fine adjustments');
tic
FINE_N = 4*HOP;

range = 3;

subsamples = 16;

% debug_offset = 3514;


tuned_delays = summed_delays(m+1:end);

% some_delays = summed_delays(m+1:m+20);
% tuned_delays = some_delays;

% tuned_delays = tuned_delays(debug_offset:debug_offset);


% RMSs = zeros(2*range*subsamples+1,length(tuned_delays));

% volumes = zeros(size(tuned_delays));
% volume_range = zeros(size(tuned_delays));
% A_volume = zeros(size(tuned_delays));
% B_volume = zeros(size(tuned_delays));

locations = bsxfun(@plus, [-FINE_N/2*subsamples+1:subsamples:FINE_N/2*subsamples]', [-range*subsamples:range*subsamples]);
B_upsampled = resample(B_, subsamples, 1);

for i = 1:length(tuned_delays)-m
  % i = debug_offset;

  j = 2*BIG_HOP + (i-1)*HOP;

  % d = tuned_delays(1);
  d = tuned_delays(i);

  d_int = floor(d);

  d_rem = floor(d*subsamples)/subsamples - d_int;

  left = j+1-FINE_N/2;
  right = j+FINE_N/2;

  frame_A = A_(left:right);

  % frame_B = B_(d_int+left-range:d_int+right+range);

  % i have a feeling this constant resampling is inefficient
  % B_upsampled = resample(frame_B, subsamples, 1);

  % basically a localized cross-correlation
  % sums = rms(bsxfun(@minus, frame_A, B_upsampled(subsamples*(d_int+d_rem+(left-1)+FINE_N/2)+1 +locations)));
  sums = sum(abs(bsxfun(@minus, frame_A, B_upsampled(subsamples*(d_int+d_rem+(left-1)+FINE_N/2)+1 +locations))),1);

  sums = 10*log10(sums);

  % figure(1);
  % subplot(2,1,1);
  % imagesc(Bs');
  % subplot(2,1,2);
  % imagesc(bsxfun(@minus, frame_A, Bs)');

  % [peak peak_i] = min(sums(2:end-1));
  % peak_i = peak_i + 1;
  [peaks peaks_i] = findpeaks(-sums);
  if (length(peaks) > 0)
    [peak max_peaks_i] = max(peaks);
    peak_i = peaks_i(max_peaks_i);
    if (peak_i > 1 & peak_i < length(sums))
      y = sums(peak_i-1:peak_i+1);
      dx = 0.5 * (y(1) - y(3)) / (y(1) - 2 * y(2) + y(3));
    else
      peak = sums(range*subsamples+1);
      dx = 0;
    end
  else
    peak = sums(range*subsamples+1);
    dx = 0;
  end

  % figure(2);

  % plot(sums, '-o', 'color','b');
  % hold on;
  % line([range*subsamples+1 range*subsamples+1],[max(sums) min(sums)], 'color','r')
  % line([(range-d_rem)*subsamples+1 (range-d_rem)*subsamples+1],[max(sums) min(sums)], 'color','g')
  % line([peak_i+dx peak_i+dx],[max(sums) min(sums)], 'color','b')

  % % f=fit([peak_i-1:peak_i+1]',sums([peak_i-1:peak_i+1])','poly2')
  % % plot(f)
  % hold off;

  % drawnow;

  d_fine = (d - range) + (peak_i + dx - 1)/subsamples;

  tuned_delays(i) = d_fine;
  % volumes(i) = peak;

  % volume_range(i) = max(sums) - min(sums);
  % A_volume(i) = 10*log10(rms(frame_A));
  % B_volume(i) = 10*log10(rms(frame_B));

end

% plot(summed_delays(m+1:end));
% hold on;
% plot(tuned_delays);
% hold off;
toc
%%%%%%%%% render
% I DON'T KNOW WHYYYYYY I HAVE TO PAD THE ZEROS
% max_delay = ceil(max(summed_delays))+8;
% summed_delays_smooth = resample(max_delay - summed_delays(m-1:end), HOP, 1);
tic
display('processing delays')

max_delay = ceil(max(tuned_delays))+8;
summed_delays_smooth = resample(max_delay - [0; 0; tuned_delays(1:end)], HOP, 1);

A_delayed = [zeros(max_delay,1); A_(2*BIG_HOP+1:end-2*BIG_HOP)];
B_trimmed = B_(2*BIG_HOP+1:end-2*BIG_HOP);

hsr = dsp.SignalSource(B_trimmed, HOP);
hms = dsp.SignalSource(summed_delays_smooth, HOP);
% hms = dsp.SignalSource(summed_delays_smooth + 0.5*sin(2*pi*[0:length(summed_delays_smooth)-1]'/(5*fs)), HOP);
% hms = dsp.SignalSource(summed_delays_smooth, HOP);
hvfd = dsp.VariableFractionalDelay('MaximumDelay', max_delay, 'InterpolationMethod', 'FIR');
hLog = dsp.SignalSink;


for i = 1:length(tuned_delays)
% for i = 1:1000
  step(hLog, step(hvfd, step(hsr), step(hms)));
end

C = hLog.Buffer;
toc
% offset = 905; % idk man

player = audioplayer( A_delayed - C(1:length(A_delayed)) , fs);
% player = audioplayer( A_delayed(1:length(C)) - C , fs);

% summed_output = A_delayed - C(1:length(A_delayed));
% rms(summed_output)

tuned_output = A_delayed - C(1:length(A_delayed));
rms(tuned_output)

play(player);
display('done');
toc
audiowrite('./audio/output.wav',tuned_output/max(max(tuned_output),-min(tuned_output)), fs);
