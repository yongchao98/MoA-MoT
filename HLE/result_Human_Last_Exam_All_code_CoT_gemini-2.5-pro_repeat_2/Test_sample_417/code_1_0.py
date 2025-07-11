# Plan:
# 1. Define the number of signals from the central core methyl groups.
# 2. Define the number of signals from the methylene linkers.
# 3. Define the number of signals from one of the identical substituent arms by summing its components.
# 4. Sum these counts to get the total number of expected 1H NMR peaks.
# 5. Print the detailed breakdown and the final calculation.

# 1. Signals from the central 2,4,6-trimethylbenzene core
# The three methyl groups are equivalent due to C3 symmetry.
core_methyl_signals = 1

# 2. Signals from the -CH2- linkers
# The three -CH2- groups are equivalent due to C3 symmetry.
# However, each CH2 is attached to a chiral center, making its two protons diastereotopic.
linker_ch2_signals = 2

# 3. Signals from one substituent arm ((4S,7R)-...-indazol-2-yl)
# This rigid, chiral moiety has several unique proton environments.
pyrazole_ring_H = 1      # One proton on the pyrazole ring
bridgehead_H = 1         # One proton on the bicyclic framework's bridgehead
c5_methylene_protons = 2 # Two diastereotopic protons
c6_methylene_protons = 2 # Two different diastereotopic protons
bridgehead_methyl = 1    # One unique methyl group
gem_dimethyls = 2        # Two diastereotopic methyl groups

# Sum the signals from one arm
substituent_arm_signals = pyrazole_ring_H + bridgehead_H + c5_methylene_protons + c6_methylene_protons + bridgehead_methyl + gem_dimethyls

# 4. Calculate the total number of signals
total_signals = core_methyl_signals + linker_ch2_signals + substituent_arm_signals

# 5. Print the analysis and result
print("Calculation of Expected 1H NMR Peaks:")
print("-" * 40)
print(f"Signals from central core methyl groups: {core_methyl_signals}")
print(f"Signals from diastereotopic -CH2- linkers: {linker_ch2_signals}")
print(f"Signals from one substituent arm: {substituent_arm_signals}")
print("  - Pyrazole ring H: 1")
print("  - Bridgehead H: 1")
print("  - C5-Methylene Hs: 2")
print("  - C6-Methylene Hs: 2")
print("  - Bridgehead Methyl: 1")
print("  - Gem-Dimethyls: 2")
print("-" * 40)
print("Total Expected Peaks = (Core Methyls) + (Linkers) + (Substituent Arm)")
print(f"Total Expected Peaks = {core_methyl_signals} + {linker_ch2_signals} + {substituent_arm_signals} = {total_signals}")
