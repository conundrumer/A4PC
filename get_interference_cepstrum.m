function [ out ] = get_interference_cepstrum( A, B, frame_size, frame_hop, cepstrum_length_ratio )

% returns a cepstrum spectrogram of the interference between A and B
NOISE_FLOOR = 0.001; % -30 dB

function [spec] = get_spectrogram(x)
  spec = abs(spectrogram(x, frame_size, frame_size-frame_hop, 'yaxis'));
  spec = spec(1:frame_size/2,:); % exclude nyquist
end

function [spec_diff] = get_spec_difference(x, y)
  spec_diff = log(max(NOISE_FLOOR, x - y));
end

A_spectrogram = get_spectrogram(A);

add_diff = get_spec_difference(A_spectrogram, get_spectrogram(A + B));
sub_diff = get_spec_difference(A_spectrogram, get_spectrogram(A - B));

fft_size = cepstrum_length_ratio * frame_size/2;

out = abs(fft(add_diff - sub_diff, fft_size, 1));
out = out(1:fft_size/2, :);

end

