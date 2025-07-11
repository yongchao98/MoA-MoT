# Step 1: Define the number of signals from the central trimethylbenzene core.
# The three methyl groups at positions 2, 4, and 6 are equivalent due to C3 symmetry.
central_methyl_signals = 1

# Step 2: Define the number of signals from the three linker methylene (-CH2-) groups.
# The three -CH2- groups are equivalent by C3 symmetry.
# However, each CH2 is attached to a chiral center, making its two protons diastereotopic.
# This results in two distinct signals for the methylene protons.
linker_methylene_signals = 2

# Step 3: Define the number of signals from one of the indazole substituent units.
# This unit is rigid and chiral, leading to many non-equivalent protons.
# H-3 proton on the pyrazole ring
h3_signals = 1
# H-4 bridgehead proton
h4_signals = 1
# H-5 methylene protons (diastereotopic)
h5_signals = 2
# H-6 methylene protons (diastereotopic)
h6_signals = 2
# C7-methyl group
c7_methyl_signals = 1
# C8-gem-dimethyl groups (diastereotopic)
c8_dimethyl_signals = 2

# Sum the signals from the single substituent unit
substituent_signals = h3_signals + h4_signals + h5_signals + h6_signals + c7_methyl_signals + c8_dimethyl_signals

# Step 4: Calculate the total number of signals by summing the signals from all parts.
total_signals = central_methyl_signals + linker_methylene_signals + substituent_signals

# Step 5: Print the breakdown and the final result.
print("Calculation of 1H NMR Peaks:")
print(f"Signals from central core methyl groups: {central_methyl_signals}")
print(f"Signals from linker methylene (-CH2-) groups: {linker_methylene_signals}")
print(f"Signals from one indazole unit: {h3_signals} (H-3) + {h4_signals} (H-4) + {h5_signals} (H-5) + {h6_signals} (H-6) + {c7_methyl_signals} (C7-Me) + {c8_dimethyl_signals} (C8-diMe) = {substituent_signals}")
print("---")
print(f"Total expected peaks = {central_methyl_signals} + {linker_methylene_signals} + {substituent_signals} = {total_signals}")
