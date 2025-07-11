# Plan:
# 1. Account for the number of unique proton signals from the central symmetric core.
# 2. Account for the number of unique proton signals from one of the three identical substituent arms.
# 3. Sum these numbers to find the total number of expected peaks in the 1H NMR spectrum.

# The central core is 1,3,5-substituted-2,4,6-trimethylbenzene.
# The three methyl groups are equivalent due to C3 symmetry.
core_peaks = 1

# One substituent arm consists of a -CH2- linker and a camphor-indazole unit.
# The two protons of the -CH2- linker are diastereotopic because they are next to a chiral center.
linker_peaks = 2

# The camphor-indazole unit is a rigid, chiral system. We count its non-equivalent protons:
# - 1 proton on the pyrazole ring.
# - 3 non-equivalent methyl groups.
# - 5 non-equivalent protons on the saturated bicyclic frame (1 CH and 2x CH2, where each CH2 has diastereotopic protons).
indazole_unit_peaks = 1 + 3 + (1 + 2 + 2)

# Calculate the total number of peaks
total_peaks = core_peaks + linker_peaks + indazole_unit_peaks

print("Calculation of 1H NMR Peaks:")
print("1. Central Core Peaks:")
print(f"   - The 3 equivalent methyl groups on the benzene ring produce {core_peaks} peak.")
print("\n2. Substituent Arm Peaks (analyzing one arm due to C3 symmetry):")
print(f"   - The 2 diastereotopic protons of the -CH2- linker produce {linker_peaks} peaks.")
print(f"   - The rigid, chiral camphor-indazole unit produces {indazole_unit_peaks} peaks.")
print("\n3. Total Number of Peaks:")
print(f"   - The total number of expected peaks is the sum of these contributions.")
print(f"   - Total Peaks = (Core Methyls) + (Linker CH2) + (Camphor-Indazole Unit)")
print(f"   - Total Peaks = {core_peaks} + {linker_peaks} + {indazole_unit_peaks} = {total_peaks}")
