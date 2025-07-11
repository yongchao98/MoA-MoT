import collections

def calculate_peptide_fragment_mw():
    """
    Calculates the monoisotopic molecular weight of a peptide fragment containing
    an unnatural amino acid.
    """
    # Monoisotopic masses of amino acid residues (the mass of the amino acid
    # once it's part of a peptide chain, i.e., after losing H2O).
    residue_masses = {
        'A': 71.03711,   # Alanine
        'R': 156.10111,  # Arginine
        'C': 103.00919,  # Cysteine
        'G': 57.02146,   # Glycine
        'I': 113.08406,  # Isoleucine
        'K': 128.09496,  # Lysine
        'L': 113.08406,  # Leucine
        'T': 101.04768,  # Threonine
        'V': 99.06841,   # Valine
        # Unnatural Amino Acid: Azido Phenylalanine (p-azido-L-phenylalanine)
        # Full monoisotopic mass: 206.0800 Da
        # Mass of H2O: 18.01056 Da
        # Residue mass = 206.0800 - 18.01056 = 188.06944 Da
        'X': 188.06944,
    }
    
    water_mass = 18.01056
    
    # The known peptide fragment from the user's prompt
    fragment_seq = "KAVCLXVIGATR"
    
    total_residue_mass = 0
    mass_equation_parts = []
    
    print(f"Calculating Molecular Weight for fragment: {fragment_seq}\n")
    print("Individual Residue Masses (Da):")
    
    # Calculate the sum of residue masses
    for aa in fragment_seq:
        mass = residue_masses.get(aa, 0)
        print(f"  - {aa}: {mass}")
        total_residue_mass += mass
        mass_equation_parts.append(str(mass))

    # The total MW of a peptide is the sum of its residue masses plus the mass of one water molecule
    # (for the H- on the N-terminus and -OH on the C-terminus).
    total_mw = total_residue_mass + water_mass
    
    # Construct the final equation string
    equation_str = " + ".join(mass_equation_parts) + f" + {water_mass} (H2O)"
    
    print("\nFinal Calculation:")
    print(f"{equation_str} = {total_mw:.4f} Da")

    print(f"\nThe total monoisotopic molecular weight of the fragment '{fragment_seq}' is {total_mw:.4f} Da.")

# Run the calculation
calculate_peptide_fragment_mw()