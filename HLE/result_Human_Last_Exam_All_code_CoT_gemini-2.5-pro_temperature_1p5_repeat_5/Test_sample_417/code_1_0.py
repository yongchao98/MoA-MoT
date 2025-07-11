# Plan:
# 1. Determine the number of unique proton signals from the central mesitylene core.
# 2. Determine the number of unique proton signals from one of the three identical substituent arms.
# 3. The total number of peaks is the sum of signals from the core and one arm due to the molecule's C3 symmetry.

# 1. Signals from the central 2,4,6-trimethylbenzene core.
# The three methyl groups are equivalent by symmetry.
signals_from_core = 1
print(f"Number of signals from the central core: {signals_from_core}")

# 2. Signals from one substituent arm.
# The arm is -CH2-[chiral_group].
# The -CH2- linker protons are diastereotopic.
signals_from_linker_CH2 = 2

# The chiral group is (4S,7R)-7,8,8-trimethyl-4,5,6,7-tetrahydro-4,7-methano-2H-indazol-2-yl.
# This is a rigid, chiral system.
# Signal from the proton on the indazole ring.
signals_from_indazole_H = 1
# Signal from the C7-methyl group.
signals_from_C7_methyl = 1
# Signals from the diastereotopic C8 gem-dimethyl groups.
signals_from_C8_gem_dimethyls = 2
# Signal from the C4 bridgehead proton.
signals_from_C4_methine = 1
# Signals from the diastereotopic C5 methylene protons.
signals_from_C5_methylene = 2
# Signals from the diastereotopic C6 methylene protons.
signals_from_C6_methylene = 2

# Sum of signals from one arm.
signals_from_one_arm = (signals_from_linker_CH2 + 
                       signals_from_indazole_H + 
                       signals_from_C7_methyl + 
                       signals_from_C8_gem_dimethyls + 
                       signals_from_C4_methine + 
                       signals_from_C5_methylene + 
                       signals_from_C6_methylene)

print(f"Number of signals from one substituent arm: "
      f"{signals_from_linker_CH2} + {signals_from_indazole_H} + {signals_from_C7_methyl} + "
      f"{signals_from_C8_gem_dimethyls} + {signals_from_C4_methine} + "
      f"{signals_from_C5_methylene} + {signals_from_C6_methylene} = {signals_from_one_arm}")

# 3. Total number of signals.
total_signals = signals_from_core + signals_from_one_arm
print("\nFinal Calculation:")
print(f"Total Signals = (Signals from Core) + (Signals from one Arm)")
print(f"Total Signals = {signals_from_core} + {signals_from_one_arm} = {total_signals}")
