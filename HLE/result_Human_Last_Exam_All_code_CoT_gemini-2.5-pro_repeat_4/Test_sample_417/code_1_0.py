# Step 1: Define the number of signals from the methyl groups on the central benzene ring.
# Due to C3 symmetry, the three methyl groups are equivalent.
signals_from_benzene_methyls = 1

# Step 2: Define the number of signals from the linker -CH2- groups.
# The three -CH2- groups are equivalent due to C3 symmetry.
# However, each -CH2- is attached to a chiral substituent, making its two protons diastereotopic.
signals_from_linker_ch2 = 2

# Step 3: Define the number of signals from one camphor-derived substituent (R).
# This substituent is chiral and has no internal symmetry.
# a) Signals from the three non-equivalent methyl groups.
signals_from_R_methyls = 3
# b) Signal from the single proton on the pyrazole ring.
signals_from_R_pyrazole_H = 1
# c) Signal from the single proton on the bridgehead carbon.
signals_from_R_bridgehead_H = 1
# d) Signals from the two CH2 groups in the bicyclic system. Each has 2 diastereotopic protons.
signals_from_R_ring_CH2s = 4
# Total signals from one R group.
signals_from_R_group = signals_from_R_methyls + signals_from_R_pyrazole_H + signals_from_R_bridgehead_H + signals_from_R_ring_CH2s

# Step 4: Calculate the total number of signals by summing the contributions.
total_signals = signals_from_benzene_methyls + signals_from_linker_ch2 + signals_from_R_group

# Final output
# The final print statement shows how the total is derived from its components.
print("Calculation of expected 1H NMR peaks:")
print(f"Signals from central ring methyl groups: {signals_from_benzene_methyls}")
print(f"Signals from linker -CH2- groups: {signals_from_linker_ch2}")
print(f"Signals from one chiral substituent (R): {signals_from_R_group}")
print(f"Total expected peaks = {signals_from_benzene_methyls} + {signals_from_linker_ch2} + {signals_from_R_group} = {total_signals}")
