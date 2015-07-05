function [ A_, B_, delays_ ] = normalize_delay( A, B, N, HOP )
% normalize_delay: returns A and B approximately in phase, with zero
% padding. please make sure A and B are within 371.5 ms of each other
% 371.5 ms
RATIO = 1;
BEGINNING_OFFSET = 1;
ENDING_OFFSET = 6;

%% get_avg_and_dev: returns average and deviation of delays
function [avg, dev, delays] = get_avg_and_dev(A, B)
  cepstrum = get_interference_cepstrum(A, B, N, HOP, RATIO);
  [M, delays] = max(cepstrum,[],1);
  delays = 2*delays;
  delays = delays - 1;
  delays(1:BEGINNING_OFFSET) = delays(BEGINNING_OFFSET+1);
  delays(length(delays)-ENDING_OFFSET:length(delays)) = delays(length(delays)-ENDING_OFFSET-1);
  avg = round((mean(delays)));
  dev = round(std(delays));

  % subplot(2,1,1);
  % hold on;
  % plot(delays,'-o');
  % subplot(2,1,2);
  % hold on;
  % plot(M);

end

[init_avg init_dev delays] = get_avg_and_dev( A, B );

% init_avg

% max(delays)

% min(delays)

range = max(max(delays) - init_avg, init_avg - min(delays));

deviation = round(2*std(diff(delays)));

shift = init_avg - range - deviation;

padding = zeros(shift, 1 );

A_padded = [ A; padding ];
B_padded = [ B; padding ];
A_delayed = [ padding; A ];
B_delayed = [ padding; B ];

[A_delayed_avg A_delayed_dev A_delayed_delays ] = get_avg_and_dev(A_delayed, B_padded);
[B_delayed_avg B_delayed_dev B_delayed_delays ] = get_avg_and_dev(B_delayed, A_padded);

if A_delayed_avg < init_avg
  if A_delayed_avg < B_delayed_avg
    % 'A delayed by this much'
    % init_avg
    A_ = A_delayed;
    B_ = B_padded;
    delays_ = A_delayed_delays;
  end
elseif B_delayed_avg < init_avg
  if B_delayed_avg < A_delayed_avg
    % 'B delayed by this much'
    % init_avg
    A_ = A_padded;
    B_ = B_delayed;
    delays_ = B_delayed_delays;
  end
else
  % 'A and B already approximately in phase'
  A_ = A;
  B_ = B;
  delays_ = delays;
end

end

