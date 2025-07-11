# --- Parameters ---
# Efficiency of each amino acid coupling step in SPPS
spps_step_efficiency = 0.99

# Efficiency of the Native Chemical Ligation step
ligation_efficiency = 0.70

# Length of the peptides
full_length = 100
fragment_length = 50

# --- Calculations ---

# 1. Direct SPPS synthesis of the 100aa peptide
# A peptide of length N requires N-1 coupling steps.
direct_synthesis_yield = spps_step_efficiency ** (full_length - 1)

# 2. NCL approach
# First, synthesize a 50aa fragment. This is the yield-limiting synthesis step.
fragment_synthesis_yield = spps_step_efficiency ** (fragment_length - 1)

# The overall yield is the yield of the fragment synthesis multiplied by the ligation efficiency.
# We assume the second fragment is produced in sufficient quantity.
ncl_overall_yield = fragment_synthesis_yield * ligation_efficiency

# --- Output Results ---
print("--- Yield Comparison: Direct SPPS vs. Native Chemical Ligation (NCL) ---")
print(f"Assuming an SPPS coupling efficiency of {spps_step_efficiency * 100:.1f}%\n")

print("1. Direct Synthesis (100 amino acids):")
print(f"   Equation: Yield = Coupling_Efficiency ^ (Peptide_Length - 1)")
print(f"   Calculation: {spps_step_efficiency} ^ ({full_length} - 1)")
print(f"   Theoretical Yield: {direct_synthesis_yield:.4f} or {direct_synthesis_yield * 100:.2f}%\n")

print("2. NCL Approach (2 x 50 amino acid fragments):")
print(f"   Yield of one 50aa fragment synthesis:")
print(f"   Equation: Yield = Coupling_Efficiency ^ (Fragment_Length - 1)")
print(f"   Calculation: {spps_step_efficiency} ^ ({fragment_length} - 1)")
print(f"   Fragment Yield: {fragment_synthesis_yield:.4f} or {fragment_synthesis_yield * 100:.2f}%")
print(f"\n   Overall NCL yield (assuming {ligation_efficiency * 100:.0f}% ligation efficiency):")
print(f"   Equation: Final_Yield = Fragment_Yield * Ligation_Yield")
print(f"   Calculation: {fragment_synthesis_yield:.4f} * {ligation_efficiency}")
print(f"   Final NCL Theoretical Yield: {ncl_overall_yield:.4f} or {ncl_overall_yield * 100:.2f}%\n")

print("--- Conclusion ---")
print(f"The NCL approach yields approximately {ncl_overall_yield*100:.1f}%, while direct SPPS yields only {direct_synthesis_yield*100:.1f}%.")
print("This demonstrates why NCL is the superior method for synthesizing long peptides.")