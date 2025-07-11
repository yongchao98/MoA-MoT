import sys

def calculate_peptide_mw():
    """
    Calculates the monoisotopic molecular weight of a peptide.
    The plan is to:
    1. Define a representative 100-amino acid peptide sequence containing 'X' (Azido Phenylalanine).
    2. Create a dictionary of monoisotopic molecular weights for the residues of all standard amino acids and Azido Phenylalanine.
    3. Iterate through the peptide sequence, summing the masses of each residue.
    4. Add the mass of a water molecule (H2O) to account for the N- and C-termini.
    5. Print the contribution of each amino acid to the total mass, as requested.
    6. Print the final total molecular weight.
    """
    # 1. Define the peptide sequence and the symbol for Azido Phenylalanine.
    # The sequence is a plausible 100aa example based on the user's snippet "M...[KAVCLXVIGATR[...]A".
    unnatural_aa_symbol = 'X'
    peptide_sequence = "M" + "G" * 30 + "KAVCLXVIGATR" + "S" * 56 + "A"
    
    # 2. Monoisotopic masses of amino acid RESIDUES (mass - H2O).
    residue_masses = {
        'A': 71.03711,  'R': 156.10111, 'N': 114.04293, 'D': 115.02694,
        'C': 103.00919, 'E': 129.04259, 'Q': 128.05858, 'G': 57.02146,
        'H': 137.05891, 'I': 113.08406, 'L': 113.08406, 'K': 128.09496,
        'M': 131.04049, 'F': 147.06841, 'P': 97.05276,  'S': 87.03203,
        'T': 101.04768, 'W': 186.07931, 'Y': 163.06333, 'V': 99.06841,
        unnatural_aa_symbol: 188.06981 # Monoisotopic mass of Azido Phenylalanine residue (C9H8N4O)
    }
    
    # Mass of a water molecule to add back for the terminal ends.
    water_mass = 18.01056
    
    total_mw = 0
    
    # 3. & 5. Iterate, sum, and print each step.
    print(f"Calculating molecular weight for the {len(peptide_sequence)}aa peptide:")
    print(f"Sequence: {peptide_sequence}\n")
    print("Starting equation...")
    
    # Print equation starting with mass of H2O for termini
    sys.stdout.write(f"{water_mass:.4f} (H2O)")
    total_mw += water_mass
    
    # Add mass of each residue
    for aa in peptide_sequence:
        if aa in residue_masses:
            mass = residue_masses[aa]
            total_mw += mass
            # Use sys.stdout.write to prevent newlines and build the equation on one line
            sys.stdout.write(f" + {mass:.4f} ({aa})")
        else:
            print(f"\nError: Amino acid '{aa}' not found in mass dictionary.")
            return

    # 4. & 6. Print the final result
    print(f"\n\n---------------------------------")
    print(f"Final Molecular Weight of the peptide = {total_mw:.4f} Da")

# Execute the function
calculate_peptide_mw()