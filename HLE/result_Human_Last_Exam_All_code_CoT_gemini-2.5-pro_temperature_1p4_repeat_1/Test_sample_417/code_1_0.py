# Plan:
# 1. The molecule has C3 symmetry. We count the signals for one repeating unit.
# 2. A repeating unit consists of a core methyl group, a CH2 linker, and a complex chiral arm.
# 3. Sum the signals from each part to get the total.

# Number of signals from the three equivalent methyl groups on the central benzene ring.
# Due to C3 symmetry, they are all equivalent.
signals_core_methyls = 1

# Number of signals from the protons of one CH2 linker group.
# The three CH2 linkers are equivalent by C3 symmetry.
# The two protons on a single CH2 group are diastereotopic because they are adjacent to a chiral group.
# Diastereotopic protons are non-equivalent.
signals_linker_CH2 = 2

# Number of signals from one of the three equivalent chiral arms.
# The arm is a rigid, asymmetric bicyclic system.
# We count the non-equivalent protons within one arm:
# - 1 signal from the proton on the pyrazole ring.
# - 1 signal from the bridgehead proton.
# - 4 signals from the four non-equivalent protons of the two CH2 groups in the bicyclic frame.
# - 3 signals from the three non-equivalent methyl groups.
# Total = 1 + 1 + 4 + 3 = 9
signals_arm = 9

# Calculate the total number of expected 1H NMR peaks
total_signals = signals_core_methyls + signals_linker_CH2 + signals_arm

# Print the final calculation, showing each component
print(f"Calculation of total 1H NMR signals:")
print(f"Signals from central methyl groups: {signals_core_methyls}")
print(f"Signals from one CH2 linker: {signals_linker_CH2}")
print(f"Signals from one substituent arm: {signals_arm}")
print(f"Total number of signals = {signals_core_methyls} + {signals_linker_CH2} + {signals_arm} = {total_signals}")