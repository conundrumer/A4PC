

[A, fs] = audioread('./audio/full_A.wav');
A = A(:,1) + A(:,2);

[B, fs] = audioread('./audio/full_B.wav');
B = B(:,1) + B(:,2);

% B = B(850:length(B));

A = A(1:length(B));


% offset = 1 + fs*10;
% dur = fs*30;
% A = A(offset:offset+dur);
% B = B(offset:offset+dur);

% B = -B; % oops i already inverted the B samples


[A, B] = normalize_delay(A, B);

A = padarray(A, n);
B = padarray(B, n);



% 371.5 ms
N = 2048;
% N = 256*4;
HOP = N/2;
% ratio = 16;

delays = get_int_delays(A, B, N, HOP);

% plot(delays);

C_left = zeros(length(A),1);
C_right = zeros(length(A),1);

% C_left = A;
% C_right = A;
w = hann(2*hop);
% o = -hop;
o=0;

for i = 2:length(delays)-1
  x = i * hop;
  j = delays(i);
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
% player = audioplayer(A-C_left, fs);
player = audioplayer(A-C_right, fs);

play(player);
