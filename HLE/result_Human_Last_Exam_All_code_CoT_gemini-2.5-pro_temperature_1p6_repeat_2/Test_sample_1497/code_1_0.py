import sys

def solve_violin_puzzle():
  """
  This function solves the logic puzzle about violin acoustics.
  The reasoning is as follows:

  1.  (sul ponticello): Bowing near the bridge excites more high-frequency harmonics.
      The parameter `ν` in `exp(-n/ν)` controls the strength of harmonics `n`. A larger
      `ν` means stronger high harmonics. This corresponds to a change in group `ii`.

  2.  (with a bridge mute): A mute adds mass to the bridge, damping its vibrations.
      The bridge is responsible for critical high-frequency resonances. Group `iv`
      (`μ, a_2, f_2`) models a high-frequency resonance. Adding mass to an oscillator
      lowers its resonant frequency. Therefore, this changes group `iv` and `f_2`
      goes down.

  3.  (in a room filled with helium): The air inside the violin body acts as a
      resonator. The speed of sound in helium is much higher than in air, which
      will significantly increase the frequency of the main air resonance. This
      low-frequency body resonance is modeled by group `iii` (`a_1, f_1`).

  4.  (on the E string): The standard setup is on the A string. Changing to the E
      string means changing the fundamental frequency of the open string itself. This
      is represented by the parameter `F`, which is group `i`.

  Combining these findings in the specified order (1, 2, 3, 4) and adding the
  direction of change for parameter `f_2` from variation (2) gives the final answer.
  """
  # The matched groups for variations (1), (2), (3), (4)
  group_1 = "ii"
  group_2 = "iv"
  group_3 = "iii"
  group_4 = "i"

  # The direction of change for the last parameter (f2) in group_2
  direction = "down"

  # Construct the final answer string
  answer = f"{group_1},{group_2},{group_3},{group_4},{direction}"
  print(answer)

solve_violin_puzzle()