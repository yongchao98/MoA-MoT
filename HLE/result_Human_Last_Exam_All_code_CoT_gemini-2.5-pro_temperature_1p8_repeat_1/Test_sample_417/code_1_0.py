# This script calculates the number of expected 1H NMR signals for the given molecule.
# The method involves analyzing the molecule's symmetry and counting the number of
# chemically non-equivalent proton environments.

# 1. Define signals from the central core of the molecule.
# The molecule has C3 symmetry, making the 3 methyl groups on the benzene ring equivalent.
signals_from_core_methyls = 1

# The 3 -CH2- linker groups are also equivalent. However, each -CH2- group is attached
# to a chiral center, making its two protons diastereotopic. This results in two distinct signals.
signals_from_linkers = 2

# 2. Define signals from one of the three identical substituent arms.
# The arm is a rigid, chiral structure without internal symmetry. We count each unique proton environment.
#   - 1 signal from the proton on the pyrazole ring.
#   - 1 signal from the bridgehead proton.
#   - 4 signals from the four non-equivalent protons of the two CH2 groups in the saturated ring.
#   - 1 signal from the methyl group on the bridgehead.
#   - 2 signals from the two diastereotopic geminal methyl groups.
# The total number of signals from one arm is the sum of these individual signals.
signals_from_one_arm = 1 + 1 + 4 + 1 + 2

# 3. Calculate the total number of signals by summing the contributions from each part.
total_signals = signals_from_core_methyls + signals_from_linkers + signals_from_one_arm

# 4. Print the breakdown of the calculation and the final result.
print("Calculation of total expected 1H NMR signals:")
print(f"Signals from equivalent core methyl groups: {signals_from_core_methyls}")
print(f"Signals from diastereotopic linker CH2 groups: {signals_from_linkers}")
print(f"Signals from one substituent arm: {signals_from_one_arm}")
print("---")
print(f"Total signals = {signals_from_core_methyls} (core methyls) + {signals_from_linkers} (linkers) + {signals_from_one_arm} (one arm) = {total_signals}")
