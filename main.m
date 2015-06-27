n = 2048;
% n = 256*4;
hop = 512;
% ratio = 16;
ratio = 1;

[A, fs] = audioread('./audio/wiggle_A.wav');
A = A(:,1) + A(:,2);
A = padarray(A, n);

[B, fs] = audioread('./audio/wiggle_B.wav');
B = B(:,1) + B(:,2);
B = padarray(B, n);
B = -B; % oops i already inverted the B samples

cepstrum = get_interference_cepstrum(A, B, n, hop, ratio);

% plot(sum(cepstrum,2));

[M, I] = max(cepstrum,[],1);
middle = median(I);
((middle/ratio) - 1) * 2
range = 5 * ratio;
low = max(1, middle - range);
high = min(ratio*n/2, middle + range);

subplot(3,1,1);
imagesc(cepstrum(low:high,:));
subplot(3,1,2);
plot((I-1)/ratio);
subplot(3,1,3);
plot(M);

C_left = zeros(length(A),1);
C_right = zeros(length(A),1);

% C_left = A;
% C_right = A;
w = hann(2*hop);
% o = -hop;
o=0;

for i = 2:length(I)-1
  x = i * hop;
  j = 2*round((I(i) - 1)/ratio);
  if o+x-j <= 0
    % I(i) = I(i+1);
    continue;
  end
  if o+x+j > length(A)-2*hop
    % I(i) = I(i-1);
    continue;
  end
  % j = I(i) - 1;
  % if ((o+x+j-2*hop+1) < 0) or ((o+x+j) > (length(A)-2*hop))
  %   continue;
  % end

  % size(B(x+j-hop+1:x+j))
  % size(w)
  % size(C_left(x:x+hop))

  C_right(o+x:o+x+2*hop-1) = C_right(o+x:o+x+2*hop-1) + B(o+x+j:o+x+j+2*hop-1) .* w;
  j = -j;
  C_left(o+x:o+x+2*hop-1) = C_left(o+x:o+x+2*hop-1) + B(o+x+j:o+x+j+2*hop-1) .* w;
end

% soundsc(A-B, fs);
soundsc(A-C_left, fs);
% soundsc(A-C_right, fs);
