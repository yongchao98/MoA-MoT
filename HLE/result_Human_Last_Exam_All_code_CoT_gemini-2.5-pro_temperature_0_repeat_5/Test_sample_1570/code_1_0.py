import math

def calculate_peptide_mass(sequence):
    """Calculates the monoisotopic mass of a peptide sequence."""
    residue_mass = {
        'A': 71.03711, 'R': 156.10111, 'N': 114.04293, 'D': 115.02694,
        'C': 103.00919, 'E': 129.04259, 'Q': 128.05858, 'G': 57.02146,
        'H': 137.05891, 'I': 113.08406, 'L': 113.08406, 'K': 128.09496,
        'M': 131.04049, 'F': 147.06841, 'P': 97.05276, 'S': 87.03203,
        'T': 101.04768, 'W': 186.07931, 'Y': 163.06333, 'V': 99.06841
    }
    H2O_MASS = 18.01056
    
    mass = sum(residue_mass[aa] for aa in sequence)
    return round(mass + H2O_MASS, 3)

# Define constants
H_MASS = 1.007825
# Mass change for deamidation of N->D or Q->E is approx +0.984 Da
DEAMIDATION_MASS_CHANGE = 0.984

# Peptides for the second disulfide bridge, assuming one missed cleavage
peptide_1_seq = "FLIPNACSQAESK"
peptide_2_seq = "ACSLAKTAFDEAIAELDTLSEESYK"

# Calculate individual peptide masses
mass_peptide_1 = calculate_peptide_mass(peptide_1_seq)
mass_peptide_2 = calculate_peptide_mass(peptide_2_seq)

# Account for two deamidations (one N and one Q in peptide 1)
total_deamidation_mass = round(2 * DEAMIDATION_MASS_CHANGE, 3)

# Calculate the neutral mass of the cross-linked, modified peptide
# Mass = Mass(Pep1) + Mass(Pep2) + Mass(Mod) - 2*H (for disulfide bond)
linked_peptide_neutral_mass = round(mass_peptide_1 + mass_peptide_2 + total_deamidation_mass - (2 * H_MASS), 3)

# Calculate the m/z for a plausible charge state (z=4)
# z=4 is plausible due to 4 basic residues (3xK, 1xN-terminus) at low pH
charge = 4
mz_value = round((linked_peptide_neutral_mass + (charge * H_MASS)) / charge, 3)

# Print the final equation and result
print("This calculation is for the second disulfide bridge with one missed cleavage and two deamidations.")
print(f"Peptide 1: {peptide_1_seq}")
print(f"Peptide 2: {peptide_2_seq}")
print("\nEquation for the m/z value:")
print(f"m/z = (Mass(Peptide 1) + Mass(Peptide 2) + Mass(Deamidations) - 2*Mass(H) + z*Mass(H)) / z")
print(f"m/z = ({mass_peptide_1} + {mass_peptide_2} + {total_deamidation_mass} - 2*{round(H_MASS, 3)} + {charge}*{round(H_MASS, 3)}) / {charge}")
print(f"m/z = ({linked_peptide_neutral_mass} + {round(charge * H_MASS, 3)}) / {charge}")
print(f"m/z = {mz_value}")
