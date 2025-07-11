# Step-by-step calculation of the number of 1H NMR signals.

# 1. Signals from the central 2,4,6-trimethylbenzene core.
# The three methyl groups are equivalent by C3 symmetry.
core_methyl_signals = 1

# 2. Signals from the -CH2- linking group.
# The two protons are diastereotopic because they are adjacent to a chiral center.
linking_methylene_signals = 2

# 3. Signals from one chiral substituent arm.
# This substituent is derived from camphor and is a rigid chiral bicyclic system.

# 3a. The three methyl groups on the substituent are all non-equivalent.
substituent_methyl_signals = 3

# 3b. The single proton on the pyrazole part of the indazole ring.
indazole_ch_signal = 1

# 3c. The remaining non-methyl protons on the saturated bicyclic skeleton.
# This consists of one bridgehead CH (1 signal) and two CH2 groups (2 signals each).
# The protons within each CH2 group are diastereotopic.
skeleton_ch_signals = 1 + 2 + 2

# 4. Summing all the signals.
total_signals = core_methyl_signals + linking_methylene_signals + substituent_methyl_signals + indazole_ch_signal + skeleton_ch_signals

print("Calculation of 1H NMR signals:")
print(f"Signals from central core methyls: {core_methyl_signals}")
print(f"Signals from linking -CH2- group: {linking_methylene_signals}")
print(f"Signals from substituent methyls: {substituent_methyl_signals}")
print(f"Signal from indazole C-H: {indazole_ch_signal}")
print(f"Signals from substituent skeleton C-H protons: {skeleton_ch_signals}")
print("---")
print("Total number of signals = (core) + (linker) + (subst. Me) + (indazole CH) + (skeleton CHs)")
print(f"Final equation: {core_methyl_signals} + {linking_methylene_signals} + {substituent_methyl_signals} + {indazole_ch_signal} + {skeleton_ch_signals} = {total_signals}")
print("\nNumber of expected peaks:")
print(total_signals)