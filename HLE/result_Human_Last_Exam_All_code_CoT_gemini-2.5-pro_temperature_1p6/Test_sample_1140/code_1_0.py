from collections import Counter

def calculate_peptide_mw():
    """
    Calculates and prints the molecular weight of a peptide with an unnatural amino acid.
    """
    # Monoisotopic mass of the residue (the AA minus H2O).
    residue_mw = {
        'A': 71.03711,   # Alanine
        'R': 156.10111,  # Arginine
        'N': 114.04293,  # Asparagine
        'D': 115.02694,  # Aspartic acid
        'C': 103.00919,  # Cysteine
        'E': 129.04259,  # Glutamic acid
        'Q': 128.05858,  # Glutamine
        'G': 57.02146,   # Glycine
        'H': 137.05891,  # Histidine
        'I': 113.08406,  # Isoleucine
        'L': 113.08406,  # Leucine
        'K': 128.09496,  # Lysine
        'M': 131.04049,  # Methionine
        'F': 147.06841,  # Phenylalanine
        'P': 97.05276,   # Proline
        'S': 87.03203,   # Serine
        'T': 101.04768,  # Threonine
        'W': 186.07931,  # Tryptophan
        'Y': 163.06333,  # Tyrosine
        'V': 99.06841,   # Valine
        # Unnatural amino acid: p-Azido-L-phenylalanine (AzF)
        # Residue mass = Monoisotopic MW(AzF) - Monoisotopic MW(H2O)
        # 206.0798 - 18.01056 = 188.06924 Da
        'X': 188.06924
    }
    
    # Representative 100 amino acid sequence containing M...KAVCLXVIGATR...A
    peptide_sequence = (
        "MSDEKRFNTILWHGRYQCMKW"  # 21
        "KAVCLXVIGATR"           # 12 (with UAA 'X')
        "PWETGNHVIYLSPADNVQAGF"  # 21
        "KTLYICGGHERWQPGMLDHIN"  # 21
        "TSEWCNAVFPQKVWSTGYMVA"  # 25
    )

    # The weight of the termini (H- and -OH), equivalent to one water molecule.
    water_mw = 18.01056

    print(f"Calculating monoisotopic molecular weight for the {len(peptide_sequence)}aa peptide:")
    print(f"{peptide_sequence}\n")
    print("The total molecular weight is the sum of all residue weights plus the weight of one water molecule for the termini.\n")
    print("--- Calculation Breakdown ---")

    # Count the occurrences of each amino acid in the sequence
    aa_counts = Counter(peptide_sequence)
    
    # Sort the items for a consistent, readable output
    sorted_counts = sorted(aa_counts.items())
    
    sum_of_residues_mw = 0.0
    for aa, count in sorted_counts:
        # Safely get the mw, default to 0 if not found (though all should be here)
        mw = residue_mw.get(aa, 0)
        contribution = count * mw
        sum_of_residues_mw += contribution
        print(f"Contribution from {aa} ({count}x): {count:>2} * {mw:<9.5f} = {contribution:9.5f} Da")

    total_mw = sum_of_residues_mw + water_mw
    
    print("\n--- Final Equation ---")
    print("Total MW = (Sum of all residue weights) + (Weight of H2O)")
    print(f"Total MW = {sum_of_residues_mw:.5f} + {water_mw:.5f}")
    print(f"Final Molecular Weight = {total_mw:.5f} Da")

if __name__ == '__main__':
    calculate_peptide_mw()