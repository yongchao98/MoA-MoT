def calculate_peptide_molecular_weight():
    """
    Calculates the molecular weight of a peptide containing an unnatural amino acid.
    """
    # Average molecular weights (in Daltons, Da) of amino acid residues (after loss of H2O)
    residue_weights = {
        'A': 71.08, 'R': 156.19, 'N': 114.11, 'D': 115.09, 'C': 103.15,
        'E': 129.12, 'Q': 128.13, 'G': 57.05, 'H': 137.14, 'I': 113.16,
        'L': 113.16, 'K': 128.17, 'M': 131.20, 'F': 147.18, 'P': 97.12,
        'S': 87.08, 'T': 101.11, 'W': 186.21, 'Y': 163.18, 'V': 99.13,
        # X represents the unnatural amino acid, p-azido-L-phenylalanine (AzF)
        # Formula: C9H10N4O2, MW: 206.21 Da. Residue MW = 206.21 - 18.015 = 188.19 Da
        'X': 188.19
    }

    # A representative 100aa sequence based on the user's snippet "M...[KAVCLXVIGATR...]A"
    # Length check: 1 (M) + 43 (A) + 12 (snippet) + 43 (A) + 1 (A) = 100 aa
    peptide_sequence = 'M' + 'A' * 43 + 'KAVCLXVIGATR' + 'A' * 43 + 'A'
    
    print(f"Analyzing representative 100aa sequence:\n{peptide_sequence}\n")

    total_residue_weight = 0.0
    aa_count = 0

    # Sum the weights of all residues in the sequence
    for aa in peptide_sequence:
        if aa in residue_weights:
            total_residue_weight += residue_weights[aa]
            aa_count += 1
        else:
            print(f"Warning: Amino acid '{aa}' not found in weight dictionary. Skipping.")
            continue

    # Add the weight of a water molecule (H2O, ~18.015 Da) to account for the
    # free N-terminus (H) and C-terminus (OH) of the full peptide.
    water_weight = 18.015
    final_mw = total_residue_weight + water_weight

    print("--- Molecular Weight Calculation ---")
    print(f"Number of amino acids: {aa_count}")
    print(f"Sum of all residue weights: {total_residue_weight:.3f} Da")
    print(f"Weight of terminal H2O: + {water_weight:.3f} Da")
    print("------------------------------------")
    # Final equation showing each number
    print(f"Final Equation: {total_residue_weight:.3f} Da + {water_weight:.3f} Da = {final_mw:.3f} Da")
    print(f"Predicted Average Molecular Weight: {final_mw:.3f} Da")

calculate_peptide_molecular_weight()