function [ A_, B_ ] = normalize_delay( A, B )
% normalize_delay: returns A and B approximately in phase, with zero
% padding. please make sure A and B are within 371.5 ms of each other
% 371.5 ms
N = 2048*8;
HOP = 8*N;
RATIO = 1;

%% get_avg_and_dev: returns average and deviation of delays
function [avg, dev] = get_avg_and_dev(A, B)
  cepstrum = get_interference_cepstrum(A, B, N, HOP, RATIO);
  [M, delays] = max(cepstrum,[],1);
  avg = round((mean(delays) - 1) * 2);
  dev = std(delays);
end

[init_avg init_dev] = get_avg_and_dev( A, B );

padding = zeros(init_avg, 1 );

A_padded = [ A; padding ];
B_padded = [ B; padding ];
A_delayed = [ padding; A ];
B_delayed = [ padding; B ];

[A_delayed_avg A_delayed_dev] = get_avg_and_dev(A_delayed, B_padded);
[B_delayed_avg B_delayed_dev] = get_avg_and_dev(B_delayed, A_padded);

if A_delayed_avg < B_delayed_avg & A_delayed_dev < 3 * init_dev
  % 'A delayed by this much'
  % init_avg
  A_ = A_delayed;
  B_ = B_padded;
else
  % 'B delayed by this much'
  % init_avg
  A_ = A_padded;
  B_ = B_delayed;
end

end

