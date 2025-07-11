import sys

# Define monoisotopic masses for amino acid residues, water, H, and proton.
# The problem mentions rounding to the third decimal place during calculation.
# We will use precise values and round the final result.
residue_masses = {
    'A': 71.03711, 'R': 156.10111, 'N': 114.04293, 'D': 115.02694,
    'C': 103.00919, 'E': 129.04259, 'Q': 128.05858, 'G': 57.02146,
    'H': 137.05891, 'I': 113.08406, 'L': 113.08406, 'K': 128.09496,
    'M': 131.04049, 'F': 147.06841, 'P': 97.05276, 'S': 87.03203,
    'T': 101.04768, 'W': 186.07931, 'Y': 163.06333, 'V': 99.06841
}
H2O_MASS = 18.01056
H_ATOM_MASS = 1.007825
PROTON_MASS = 1.007276

def calculate_peptide_mass(sequence):
    """Calculates the neutral monoisotopic mass of a peptide sequence."""
    mass = sum(residue_masses[aa] for aa in sequence)
    return mass + H2O_MASS

def calculate_disulfide_linked_mz(pep1_seq, pep2_seq, charge):
    """Calculates the m/z of a disulfide-linked peptide pair for a given charge."""
    mass_pep1 = calculate_peptide_mass(pep1_seq)
    mass_pep2 = calculate_peptide_mass(pep2_seq)
    
    # Mass of linked peptides = sum of individual peptide masses - mass of 2 H atoms
    linked_mass_neutral = mass_pep1 + mass_pep2 - (2 * H_ATOM_MASS)
    
    # Calculate m/z = (M_neutral + z*H+)/z
    mz_value = (linked_mass_neutral + charge * PROTON_MASS) / charge
    
    return mz_value, mass_pep1, mass_pep2, linked_mass_neutral

# Peptides for the two disulfide bridges
bridge1_pep1_seq = "YDDMAACMK"
bridge1_pep2_seq = "TQGCDEAEAGEGGEN"

bridge2_pep1_seq = "FLIPNACSQAESK"
bridge2_pep2_seq = "ACSLAK"

# --- Bridge 1 Calculation ---
print("--- Analysis of Disulfide Bridge 1 ---")
print(f"Peptide 1: {bridge1_pep1_seq}")
print(f"Peptide 2: {bridge1_pep2_seq}")
# Calculate for common charge states
for z in [2, 3, 4]:
    mz, m1, m2, lm = calculate_disulfide_linked_mz(bridge1_pep1_seq, bridge1_pep2_seq, z)
    if z == 2: # Print detailed mass breakdown once
        print(f"Mass of Peptide 1: {m1:.3f}")
        print(f"Mass of Peptide 2: {m2:.3f}")
        print(f"Neutral mass of linked peptides: {lm:.3f}")
    print(f"Expected m/z for charge +{z}: {mz:.3f}")

print("\n" + "="*40 + "\n")

# --- Bridge 2 Calculation ---
print("--- Analysis of Disulfide Bridge 2 ---")
print(f"Peptide 1: {bridge2_pep1_seq}")
print(f"Peptide 2: {bridge2_pep2_seq}")
# Calculate for common charge states
for z in [2, 3, 4]:
    mz, m1, m2, lm = calculate_disulfide_linked_mz(bridge2_pep1_seq, bridge2_pep2_seq, z)
    if z == 2: # Print detailed mass breakdown once
        print(f"Mass of Peptide 1: {m1:.3f}")
        print(f"Mass of Peptide 2: {m2:.3f}")
        print(f"Neutral mass of linked peptides: {lm:.3f}")
    print(f"Expected m/z for charge +{z}: {mz:.3f}")

print("\n--- Conclusion ---")
print("Based on standard tryptic digestion, none of the calculated m/z values for the two disulfide bridges match the provided answer choices.")
print("This suggests a possible modification, missed cleavage, or an error in the problem description.")
print("However, if forced to choose, further complex analysis involving scenarios like missed cleavages and common modifications would be required. Without a clear logical path to any of the options, any selection would be speculative.")
