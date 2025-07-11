import math

# Step 1: Define the constants for mass calculation
# Monoisotopic masses of amino acids
aa_masses = {
    'A': 71.03711, 'R': 156.10111, 'N': 114.04293, 'D': 115.02694,
    'C': 103.00919, 'E': 129.04259, 'Q': 128.05858, 'G': 57.02146,
    'H': 137.05891, 'I': 113.08406, 'L': 113.08406, 'K': 128.09496,
    'M': 131.04049, 'F': 147.06841, 'P': 97.05276, 'S': 87.03203,
    'T': 101.04768, 'W': 186.07931, 'Y': 163.06333, 'V': 99.06841
}
# Other relevant masses
mass_H2O = 18.010565      # Mass of a water molecule
mass_H = 1.007825         # Mass of a hydrogen atom
mass_proton = 1.007276    # Mass of a proton (for m/z calculation)
mass_oxidation = 15.994915 # Mass of an oxygen atom for Met oxidation

# Step 2: Define the sequences of the two peptides forming the disulfide bridge
# This scenario assumes one missed cleavage and two methionine oxidations for Bridge 1.
peptide_1 = "LAEQAERYDDMAACMK"
peptide_2_missed_cleavage = "DNLTLWTSDRTQGCDEAEAGEGGEN"

# Step 3: Calculate the mass of each individual peptide, including modifications
def calculate_peptide_mass(sequence, num_oxidations=0):
    """Calculates the neutral monoisotopic mass of a peptide."""
    residue_mass = sum(aa_masses[aa] for aa in sequence)
    mod_mass = num_oxidations * mass_oxidation
    # A peptide's mass is the sum of its residue masses plus one water molecule
    peptide_mass = residue_mass + mod_mass + mass_H2O
    return peptide_mass

# Peptide 1 has two Methionines, both assumed to be oxidized.
num_oxidations_p1 = 2
mass_p1_mod = calculate_peptide_mass(peptide_1, num_oxidations_p1)

# Peptide 2 resulted from a missed cleavage, but has no further modifications.
mass_p2_missed = calculate_peptide_mass(peptide_2_missed_cleavage)

print("Calculating mass for the first disulfide bridge with modifications...")
print(f"Peptide 1: {peptide_1}")
print(f"Number of oxidized methionines in Peptide 1: {num_oxidations_p1}")
print(f"Mass of modified Peptide 1 = (Sum of AA masses) + ({num_oxidations_p1} * Mass(O)) + Mass(H2O)")
print(f"Mass of modified Peptide 1 = {sum(aa_masses[aa] for aa in peptide_1):.3f} + ({num_oxidations_p1} * {mass_oxidation:.3f}) + {mass_H2O:.3f} = {mass_p1_mod:.3f}\n")

print(f"Peptide 2 (with one missed cleavage): {peptide_2_missed_cleavage}")
print(f"Mass of Peptide 2 = (Sum of AA masses) + Mass(H2O)")
print(f"Mass of Peptide 2 = {sum(aa_masses[aa] for aa in peptide_2_missed_cleavage):.3f} + {mass_H2O:.3f} = {mass_p2_missed:.3f}\n")

# Step 4: Calculate the mass of the disulfide-linked peptide complex
# Forming a disulfide bridge removes one Hydrogen atom from each Cysteine's thiol group.
mass_disulfide_linked = mass_p1_mod + mass_p2_missed - (2 * mass_H)
print("Calculating the neutral mass of the disulfide-linked complex...")
print("Mass_complex = Mass_peptide_1 + Mass_peptide_2 - (2 * Mass(H))")
print(f"Mass_complex = {mass_p1_mod:.3f} + {mass_p2_missed:.3f} - (2 * {mass_H:.3f}) = {mass_disulfide_linked:.3f}\n")

# Step 5: Calculate the m/z for a given charge state (z)
# For a peptide of this size, a +4 charge state is common in LC/MS.
z = 4
mz_value = (mass_disulfide_linked + (z * mass_proton)) / z

print(f"Calculating the m/z value for the complex at charge state z = +{z}...")
print("m/z = (Mass_complex + (z * Mass(proton))) / z")
print(f"m/z = ({mass_disulfide_linked:.3f} + ({z} * {mass_proton:.3f})) / {z}")
print(f"Final calculated m/z = {round(mz_value, 3)}")