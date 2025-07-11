def calculate_peptide_mass(sequence):
    """Calculates the monoisotopic mass of a peptide sequence."""
    
    # Monoisotopic residue masses of amino acids
    amino_acid_masses = {
        'A': 71.03711, 'R': 156.10111, 'N': 114.04293, 'D': 115.02694,
        'C': 103.00919, 'E': 129.04259, 'Q': 128.05858, 'G': 57.02146,
        'H': 137.05891, 'I': 113.08406, 'L': 113.08406, 'K': 128.09496,
        'M': 131.04049, 'F': 147.06841, 'P': 97.05276, 'S': 87.03203,
        'T': 101.04768, 'W': 186.07931, 'Y': 163.06333, 'V': 99.06841
    }
    
    # Mass of a water molecule to add to the sum of residue masses
    h2o_mass = 18.010565
    
    # Sum of residue masses
    mass = sum(amino_acid_masses[aa] for aa in sequence)
    
    # Add mass of terminal water
    total_mass = mass + h2o_mass
    return round(total_mass, 3)

# Define constants for calculation
hydrogen_mass = 1.007825
oxygen_mass = 15.994915
proton_mass = 1.007276

# Define peptide sequences for the hypothesized cross-link
# This involves one missed cleavage (LAEQAER + YDDMAACMK)
peptide1_sequence = "LAEQAERYDDMAACMK"
peptide2_sequence = "TQGCDEAEAGEGGEN"

# Calculate mass of each peptide
peptide1_mass = calculate_peptide_mass(peptide1_sequence)
peptide2_mass = calculate_peptide_mass(peptide2_sequence)

# Calculate the mass of the cross-linked and oxidized complex
# M = Mass(Pep1) + Mass(Pep2) - 2*Mass(H) + Mass(O)
crosslinked_mass = peptide1_mass + peptide2_mass - (2 * hydrogen_mass) + oxygen_mass
crosslinked_mass = round(crosslinked_mass, 3)


# Calculate the m/z value for a +3 charge state
charge = 3
mz_value = (crosslinked_mass + (charge * proton_mass)) / charge
mz_value = round(mz_value, 3)

# Print the detailed calculation
print("This m/z value corresponds to the first disulfide bridge with one missed cleavage and one methionine oxidation.")
print(f"Peptide 1 ('{peptide1_sequence}') mass: {peptide1_mass}")
print(f"Peptide 2 ('{peptide2_sequence}') mass: {peptide2_mass}")
print("\nCalculation of the final m/z value:")
print(f"Mass of cross-linked, oxidized peptide = Peptide 1 Mass ({peptide1_mass}) + Peptide 2 Mass ({peptide2_mass}) - 2 * H ({round(2*hydrogen_mass, 3)}) + O ({round(oxygen_mass, 3)})")
print(f"                                   = {crosslinked_mass}")
print(f"\nm/z = (Mass ({crosslinked_mass}) + {charge} * Proton Mass ({proton_mass})) / {charge}")
print(f"    = {mz_value}")