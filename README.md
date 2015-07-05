# A4PC
Automatic Analog Audio Alignment for Phase Cancellation

## Background

When you have a song, and an instrumental version of the song, and the instrumental is from the same master source as the original, theoretically, you could get the acapella by phase inverting the instrumental and mixing it with the original to phase cancel the instrumental, if you manage to align them down to the exact sample.

This works if the recordings are digital, but if they're analog, due to physical imperfection (tape stretch, wow, flutter), the speeds of the two tracks would vary just slightly enough for the two tracks to be occasionally out of phase. What you will hear, instead, is a comb filter effect varying over time, somewhat like a flanger effect (also originally done through tape delay!). You could manually stretch and warp the instrumental to be more aligned with the original, but that's a tedious and repetitive process that the computer can do by itself.

## Algorithm

1. Coarse alignment: remove large delays
2. Medium alignment: analyze delays and fit a smooth curve to describe the change in delay.
3. Fine alignment: optimize volume reduction for each frame
4. Delay modulation: Modulate the instrumental track according to the generated delays.
