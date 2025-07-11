import math

# --- Assumptions ---
peptide_length = 100
# We assume a very optimistic 99% efficiency for each coupling/deprotection cycle in SPPS.
spps_step_efficiency = 0.99
# For the NCL strategy, we split the 100aa peptide into two 50aa fragments.
fragment_1_length = 50
fragment_2_length = 50
# We assume a 70% yield for the Native Chemical Ligation step itself.
ncl_step_efficiency = 0.70

# --- Calculation for Direct SPPS of the 100aa Peptide ---
# The number of coupling steps is the length of the peptide minus 1.
num_spps_steps_direct = peptide_length - 1
# The final yield is the step efficiency raised to the power of the number of steps.
yield_direct_spps_fraction = spps_step_efficiency ** num_spps_steps_direct
yield_direct_spps_percent = yield_direct_spps_fraction * 100

print("--- Strategy 1: Direct Solid-Phase Peptide Synthesis (SPPS) ---")
print(f"Peptide Length: {peptide_length} amino acids")
print(f"SPPS Step Efficiency: {spps_step_efficiency * 100}%")
print(f"Equation: Final Yield = (Step Efficiency ^ (Length - 1)) * 100")
print(f"Numbers for Equation: ({spps_step_efficiency} ^ {num_spps_steps_direct}) * 100")
print(f"Calculated Theoretical Yield: {yield_direct_spps_percent:.2f}%\n")


# --- Calculation for Native Chemical Ligation (NCL) Strategy ---
# 1. Calculate yield for synthesizing Fragment 1 (50aa)
num_spps_steps_frag1 = fragment_1_length - 1
yield_frag1_fraction = spps_step_efficiency ** num_spps_steps_frag1

# 2. Calculate yield for synthesizing Fragment 2 (50aa)
num_spps_steps_frag2 = fragment_2_length - 1
yield_frag2_fraction = spps_step_efficiency ** num_spps_steps_frag2

# 3. Calculate final yield after ligation
# Note: For simplicity, we assume we start with molar equivalents of fragments.
# The overall yield is the product of the synthesis yields and the ligation yield.
final_yield_ncl_fraction = yield_frag1_fraction * yield_frag2_fraction * ncl_step_efficiency
final_yield_ncl_percent = final_yield_ncl_fraction * 100

print("--- Strategy 2: Native Chemical Ligation (NCL) ---")
print(f"Fragment 1 Length: {fragment_1_length} aa -> Synthesis Steps: {num_spps_steps_frag1}")
print(f"Fragment 2 Length: {fragment_2_length} aa -> Synthesis Steps: {num_spps_steps_frag2}")
print(f"NCL Step Efficiency: {ncl_step_efficiency * 100}%")
print(f"Equation: Final Yield = (Yield_Frag1 * Yield_Frag2 * NCL_Efficiency) * 100")
# To show the numbers, we'll expand the fragment yield parts
print(f"Numbers for Equation: (({spps_step_efficiency}^{num_spps_steps_frag1}) * ({spps_step_efficiency}^{num_spps_steps_frag2}) * {ncl_step_efficiency}) * 100")
print(f"Calculated Theoretical Yield: {final_yield_ncl_percent:.2f}%")
