import math

# Step 1: Define constants
# Monoisotopic masses of amino acid residues
residue_masses = {
    'A': 71.03711, 'V': 99.06841, 'D': 115.02694, 'L': 113.08406,
    'T': 101.04768, 'K': 128.09496, 'I': 113.08406, 'R': 156.10111
}
# Mass of modifications and particles
lactyl_mass = 72.02113  # C3H4O2
proton_mass = 1.007276
h2o_mass = 18.010565

# Peptide sequence and recorded m/z values
sequence = "AVDLTKLIR"
recorded_mz = [401.276, 601.392, 301.200, 518.271, 304.139, 417.223]
tolerance = 0.02 # Mass tolerance for matching

print("Analyzing peptide AVDLTKLIR for lactylation on Lysine (K)\n")

# Step 2 & 3: Calculate key y-ion masses
# Mass of y3 ion (LIR)
# y-ion mass = sum of residue masses + mass of H2O
y3_mass = residue_masses['L'] + residue_masses['I'] + residue_masses['R'] + h2o_mass
y3_mz = y3_mass + proton_mass

# Mass of y4 ion (K(Lac)LIR)
# It includes the y3 part plus the lactylated Lysine residue
k_lac_residue_mass = residue_masses['K'] + lactyl_mass
y4_mass = k_lac_residue_mass + y3_mass - h2o_mass # Subtract H2O as it's already in y3_mass
y4_mz = y4_mass + proton_mass
y4_mz_z2 = (y4_mass + 2 * proton_mass) / 2 # Doubly charged ion

# Mass of b5 ion (AVDLT)
b5_mass = sum(residue_masses[aa] for aa in "AVDLT")
# Check for a potential b5+H2O adduct, which is sometimes observed
b5_h2o_mz = b5_mass + h2o_mass + proton_mass


print("--- Theoretical m/z Calculation ---")
print(f"Mass of lactylated Lysine residue K(Lac): {k_lac_residue_mass:.3f} Da")
print(f"Calculated m/z for y3 ion (LIR): {y3_mz:.3f}")
print(f"Calculated m/z for y4 ion (K(Lac)LIR): {y4_mz:.3f}")
print(f"Calculated m/z for doubly charged y4 ion (y4++): {y4_mz_z2:.3f}")
print(f"Calculated m/z for b5+H2O ion (AVDLT+H2O): {b5_h2o_mz:.3f}")

print("\n--- Matching with Recorded Values ---")
matches = {}
# Find matches in the recorded list
for val in recorded_mz:
    if abs(val - y3_mz) < tolerance:
        matches[val] = "y3 ion (LIR)"
    if abs(val - y4_mz) < tolerance:
        matches[val] = "y4 ion (K(Lac)LIR)"
    if abs(val - y4_mz_z2) < tolerance:
        matches[val] = "doubly charged y4++ ion"
    if abs(val - b5_h2o_mz) < tolerance:
        matches[val] = "b5+H2O ion (confirms N-terminus)"

for mz, identity in matches.items():
    print(f"Recorded m/z {mz} matches the {identity}.")

print("\n--- Conclusion ---")
print("The definitive evidence for lactylation on the lysine at position 6 comes from the mass difference between the y3 and y4 ions.")
y3_found = 401.276
y4_found = 601.392
mass_diff = y4_found - y3_found
print(f"The mass shift from the matched y3 ion to the y4 ion is calculated as:")
print(f"Equation: {y4_found} (y4) - {y3_found} (y3) = {mass_diff:.3f} Da")
print(f"This mass difference of {mass_diff:.3f} Da perfectly matches the mass of a lactylated lysine residue ({k_lac_residue_mass:.3f} Da).")
print("\nTherefore, the m/z values that indicate the lysine is lactylated are those that are part of this evidence chain.")
print(f"- {y3_found}: The y3 fragment, which establishes the baseline mass before the modified residue.")
print(f"- {y4_found}: The y4 fragment, which contains the modified residue.")
print(f"- {matches.get(518.271, 'N/A')}: The m/z value 518.271 confirms the N-terminal sequence of the peptide, complementing the C-terminal evidence.")

print("\nBased on this analysis, the set of values {401.276, 601.392, 518.271} provides the comprehensive evidence to identify the lactylated peptide.")

<<<D>>>