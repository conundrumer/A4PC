[A, fs] = audioread('./audio/close_A.wav');
A = A(:,1) + A(:,2);

[B, fs] = audioread('./audio/close_B.wav');
B = B(:,1) + B(:,2);

n = 2048;
% n = 256*4;
hop = n/4;
ratio = 16;
% ratio = 1;

cepstrum = get_interference_cepstrum(A, B, n, hop, ratio);

% plot(sum(cepstrum,2));

[M, I] = max(cepstrum,[],1);
middle = median(I);
middle/ratio
range = 5 * ratio;
low = max(1, middle - range);
high = min(ratio*n/2, middle + range);

subplot(3,1,1);
imagesc(cepstrum(low:high,:));
subplot(3,1,2);
plot((I-1)/ratio);
subplot(3,1,3);
plot(M);
