import textwrap

def calculate_peptide_mw():
    """
    Calculates the molecular weight of a 100aa peptide containing azido phenylalanine.
    """
    # Residue molecular weights (mass of AA - mass of H2O).
    # The mass for azido phenylalanine (p-Azido-L-phenylalanine, AzF) is C9H10N4O2.
    # Full MW of AzF = 206.20 g/mol. Residue MW = 206.20 - 18.02 = 188.18 g/mol.
    residue_mw = {
        'A': 71.08, 'R': 156.19, 'N': 114.10, 'D': 115.09, 'C': 103.14,
        'E': 129.12, 'Q': 128.13, 'G': 57.05, 'H': 137.14, 'I': 113.16,
        'L': 113.16, 'K': 128.17, 'M': 131.19, 'F': 147.18, 'P': 97.12,
        'S': 87.08, 'T': 101.11, 'W': 186.21, 'Y': 163.18, 'V': 99.13,
        'X': 188.18  # Azido Phenylalanine (AzF)
    }
    H2O_MW = 18.02

    # Construct a plausible 100aa sequence based on the user's snippet.
    # M...[KAVCLXVIGATR[...]A
    # M + (20*A) + KAVCLXVIGATR + (66*A) + A = 100 aa
    prefix = 'A' * 20
    middle = 'KAVCLXVIGATR'
    suffix = 'A' * 66
    sequence = 'M' + prefix + middle + suffix + 'A'

    print(f"Calculating Molecular Weight for the 100aa peptide.")
    print("Assuming X = Azido Phenylalanine (AzF).")
    print("\nAssumed full sequence:")
    # Wrap the sequence text for better readability
    print('\n'.join(textwrap.wrap(sequence, width=80)))
    print("-" * 80)


    total_mw = H2O_MW
    equation_parts = [f"{H2O_MW:.2f} (H2O)"]

    for aa in sequence:
        if aa in residue_mw:
            mass = residue_mw[aa]
            total_mw += mass
            equation_parts.append(f"{mass:.2f} ({aa})")
        else:
            print(f"Warning: Amino acid '{aa}' not found in mass dictionary.")
            return

    # Print the full equation
    print("Molecular Weight Calculation Equation:")
    # The equation can be very long, so we wrap it.
    equation_str = " + ".join(equation_parts)
    print('\n'.join(textwrap.wrap(equation_str, width=80, subsequent_indent="  ")))
    
    print("-" * 80)
    print(f"Final Calculated Molecular Weight = {total_mw:.2f} Da")


if __name__ == '__main__':
    calculate_peptide_mw()