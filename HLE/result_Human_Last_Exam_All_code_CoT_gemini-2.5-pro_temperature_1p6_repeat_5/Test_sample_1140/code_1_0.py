import sys

def calculate_peptide_mw():
    """
    Calculates the molecular weight of a peptide fragment containing an unnatural amino acid.
    The final output shows the full addition equation as requested.
    """
    # Average molecular weights of amino acid residues (mass after losing H2O) in Daltons (g/mol)
    residue_mw = {
        'A': 71.0788, 'R': 156.1875, 'N': 114.1038, 'D': 115.0886,
        'C': 103.1388, 'E': 129.1155, 'Q': 128.1307, 'G': 57.0519,
        'H': 137.1411, 'I': 113.1594, 'L': 113.1594, 'K': 128.1741,
        'M': 131.1926, 'F': 147.1766, 'P': 97.1167, 'S': 87.0782,
        'T': 101.1051, 'W': 186.2132, 'Y': 163.1760, 'V': 99.1326,
    }
    # Molecular weight of p-Azido-L-phenylalanine as a residue.
    # The full molecule (C9H10N4O2) is 206.20 g/mol. As a residue (losing H2O), it's 188.185 g/mol.
    residue_mw['X'] = 188.185

    # Mass of a water molecule to be added to the sum of residue weights for the final peptide mass
    water_mw = 18.0153

    # The peptide fragment from the user prompt
    sequence = "KAVCLXVIGATR"

    total_mw = 0
    equation_parts = []

    print(f"Calculating Molecular Weight for the peptide fragment: {sequence}")
    print("-" * 50)
    print("The final mass is the sum of the residue masses plus one water molecule for the terminal groups.")
    print("Equation:")

    for aa in sequence:
        if aa in residue_mw:
            mass = residue_mw[aa]
            total_mw += mass
            equation_parts.append(f"{aa}({mass:.2f})")
        else:
            # Redirect error to stderr to not interfere with stdout if piped
            print(f"Error: Amino acid '{aa}' not found in molecular weight dictionary.", file=sys.stderr)
            return

    # Add the mass of water for the termini
    total_mw += water_mw
    equation_parts.append(f"H2O({water_mw:.2f})")

    # Construct and print the final equation and result
    equation_string = " + ".join(equation_parts)
    print(f"{equation_string} = {total_mw:.2f} g/mol")

if __name__ == '__main__':
    calculate_peptide_mw()