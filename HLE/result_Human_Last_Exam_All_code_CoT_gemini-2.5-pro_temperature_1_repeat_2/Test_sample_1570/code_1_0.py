import sys
import io

# Capture original stdout
original_stdout = sys.stdout
# Create a string buffer
string_buffer = io.StringIO()
# Redirect stdout to the string buffer
sys.stdout = string_buffer

# Monoisotopic residue masses (Da)
residue_masses = {
    'A': 71.03711, 'R': 156.10111, 'N': 114.04293, 'D': 115.02694,
    'C': 103.00919, 'E': 129.04259, 'Q': 128.05858, 'G': 57.02146,
    'H': 137.05891, 'I': 113.08406, 'L': 113.08406, 'K': 128.09496,
    'M': 131.04049, 'F': 147.06841, 'P': 97.05276, 'S': 87.03203,
    'T': 101.04768, 'W': 186.07931, 'Y': 163.06333, 'V': 99.06841
}

# Other masses (Da)
mass_H2O = 18.010565
mass_proton = 1.007276

# Peptide sequence from the first disulfide bridge locator
peptide_seq = "TQGCDEAEAGEG"

# --- Mass Calculation ---
# Sum of residue masses
sum_of_residues = 0
equation_str = []
for aa in peptide_seq:
    mass = residue_masses[aa]
    sum_of_residues += mass
    equation_str.append(f"{mass} ({aa})")

# Monoisotopic mass of the neutral peptide (M)
# M = sum of residue masses + mass of one water molecule
neutral_mass = sum_of_residues + mass_H2O

# --- m/z Calculation ---
# For a singly charged ion [M+H]+, the charge z is 1.
# m/z = (M + 1 * mass_proton) / 1
charge = 1
mz_value = (neutral_mass + charge * mass_proton) / charge

# --- Output the detailed calculation ---
print("Calculation for the m/z of peptide 'TQGCDEAEAGEG':\n")
print("1. Calculate the sum of monoisotopic residue masses:")
print(" + ".join(equation_str) + f" = {round(sum_of_residues, 3)}")

print("\n2. Calculate the neutral monoisotopic mass (M) of the peptide:")
print(f"M = (Sum of residue masses) + Mass(H2O)")
print(f"M = {round(sum_of_residues, 3)} + {mass_H2O} = {round(neutral_mass, 3)}")

print("\n3. Calculate the m/z for the singly charged ion [M+H]+:")
print("m/z = (M + Mass(proton)) / 1")
print(f"m/z = ({round(neutral_mass, 3)} + {mass_proton}) / 1 = {round(mz_value, 3)}")

# Restore original stdout
sys.stdout = original_stdout
# Get the content from the string buffer
output = string_buffer.getvalue()

# Print the captured output
print(output)