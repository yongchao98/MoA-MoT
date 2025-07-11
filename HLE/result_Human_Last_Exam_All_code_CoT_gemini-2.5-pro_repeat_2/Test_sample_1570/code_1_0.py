# Define monoisotopic masses for amino acids, water, and proton.
# Using standard values, rounded for clarity in the printout.
masses = {
    'A': 71.037, 'C': 103.009, 'D': 115.027, 'E': 129.043,
    'G': 57.021, 'Q': 128.059, 'T': 101.048
}
H2O_MASS = 18.011
PROTON_MASS = 1.007

# The peptide fragment from the problem description that leads to an answer.
peptide_sequence = "TQGCDEAEAGEG"

# Calculate the peptide mass
residue_mass_sum = 0
print("Calculating mass for peptide TQGCDEAEAGEG:")
for aa in peptide_sequence:
    mass = masses[aa]
    residue_mass_sum += mass
    print(f"Mass of {aa}: {mass}")

# The mass of a peptide is the sum of its residue masses plus one water molecule.
peptide_mass = residue_mass_sum + H2O_MASS
print(f"Sum of residue masses: {round(residue_mass_sum, 3)}")
print(f"Mass of H2O: {H2O_MASS}")
print(f"Total neutral peptide mass: {round(peptide_mass, 3)}")

# Calculate the m/z for a singly charged ion [M+H]+
charge = 1
mz_value = (peptide_mass + charge * PROTON_MASS) / charge

print("\nCalculating the m/z for the [M+H]+ ion:")
print(f"Equation: (Peptide Mass + Charge * Proton Mass) / Charge")
print(f"m/z = ({round(peptide_mass, 3)} + {charge} * {PROTON_MASS}) / {charge}")
print(f"Final m/z = {round(mz_value, 3)}")