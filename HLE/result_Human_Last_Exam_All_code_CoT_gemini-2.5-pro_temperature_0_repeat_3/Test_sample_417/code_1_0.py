# This script calculates the expected number of peaks in the 1H NMR spectrum
# of 1,3,5-tri[((4S,7R)-7,8,8-trimethyl-4,5,6,7-tetrahydro-4,7-methano-2H-indazol-2-yl)methyl]-2,4,6-trimethylbenzene
# based on its chemical structure and symmetry.

# The molecule has C3 symmetry, so we count signals for 1/3 of the molecule.

# 1. Signals from the central 2,4,6-trimethylbenzene core.
# The three methyl groups are equivalent due to C3 symmetry.
signals_ar_ch3 = 1

# 2. Signals from the -CH2- linker.
# The attached indazole group is chiral, making the two CH2 protons diastereotopic.
signals_benzylic_ch2 = 2

# 3. Signals from one (4S,7R)-7,8,8-trimethyl-4,5,6,7-tetrahydro-4,7-methano-2H-indazol-2-yl group.
# This is a rigid, chiral system, so most protons are non-equivalent.
# Proton on the indazole ring (H3).
signals_indazole_h3 = 1
# Methyl groups: one at C7, two diastereotopic at C8. Total 3 unique methyl groups.
signals_indazole_methyls = 3
# Protons on the bicyclic skeleton:
# One proton at bridgehead C4.
signals_c4_h = 1
# Two diastereotopic protons at C5.
signals_c5_ch2 = 2
# Two diastereotopic protons at C6.
signals_c6_ch2 = 2

# Calculate the total number of signals by summing the signals from each part.
total_signals = (signals_ar_ch3 +
                 signals_benzylic_ch2 +
                 signals_indazole_h3 +
                 signals_indazole_methyls +
                 signals_c4_h +
                 signals_c5_ch2 +
                 signals_c6_ch2)

# Print the breakdown of the calculation.
print("Calculation of 1H NMR signals:")
print(f"Central core methyls (x3): {signals_ar_ch3} signal")
print(f"Benzylic -CH2- linkers (x3): {signals_benzylic_ch2} signals")
print(f"Indazole ring proton (x3): {signals_indazole_h3} signal")
print(f"Indazole methyl groups (x9 total): {signals_indazole_methyls} signals")
print(f"Indazole C4-H proton (x3): {signals_c4_h} signal")
print(f"Indazole C5-H2 protons (x3): {signals_c5_ch2} signals")
print(f"Indazole C6-H2 protons (x3): {signals_c6_ch2} signals")
print("-" * 20)
print(f"Total signals = {signals_ar_ch3} + {signals_benzylic_ch2} + {signals_indazole_h3} + {signals_indazole_methyls} + {signals_c4_h} + {signals_c5_ch2} + {signals_c6_ch2}")
print(f"Total expected peaks = {total_signals}")