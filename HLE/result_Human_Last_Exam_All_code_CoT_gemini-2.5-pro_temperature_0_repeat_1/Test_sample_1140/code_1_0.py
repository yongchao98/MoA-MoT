def calculate_peptide_molecular_weight():
    """
    Calculates the monoisotopic molecular weight of a given peptide sequence.
    This example uses the fragment KAVCLXVIGATR, where X is azido phenylalanine.
    """
    # The peptide sequence fragment provided by the user
    sequence = "KAVCLXVIGATR"

    # Monoisotopic molecular weights of amino acids (in Daltons, g/mol)
    # These are the weights of the free, non-polymerized amino acids.
    aa_weights = {
        'A': 89.047679,   # Alanine
        'R': 174.111676,  # Arginine
        'C': 121.019749,  # Cysteine (assuming reduced form)
        'G': 75.032029,   # Glycine
        'I': 131.094629,  # Isoleucine
        'L': 131.094629,  # Leucine
        'K': 146.105528,  # Lysine
        'T': 119.063329,  # Threonine
        'V': 117.078979,  # Valine
        'X': 206.080376   # Azido Phenylalanine (C9H10N4O2)
    }

    # Molecular weight of a water molecule (H2O) to be subtracted for each peptide bond
    water_weight = 18.010565

    print(f"Calculating molecular weight for peptide fragment: {sequence}\n")

    # Build the equation string for display
    equation_symbols = " + ".join([f"MW({aa})" for aa in sequence])
    num_bonds = len(sequence) - 1
    print(f"Formula: ( {equation_symbols} ) - ({num_bonds} * MW(H2O))")
    print("-" * 50)

    # Build the calculation string with actual numbers
    sum_of_aa_weights = 0
    calculation_values = []
    for aa in sequence:
        weight = aa_weights.get(aa)
        if weight is None:
            print(f"Error: Amino acid '{aa}' not found in the dictionary.")
            return
        sum_of_aa_weights += weight
        calculation_values.append(f"{weight:.4f}")

    calculation_str = " + ".join(calculation_values)
    
    # Calculate the final peptide weight
    total_water_weight = num_bonds * water_weight
    final_peptide_weight = sum_of_aa_weights - total_water_weight

    # Print the final equation with all numbers
    print("Calculation:")
    print(f"( {calculation_str} ) - ({num_bonds} * {water_weight:.4f}) = Final_MW")
    print(f"{sum_of_aa_weights:.4f} - {total_water_weight:.4f} = {final_peptide_weight:.4f} Da")
    print("-" * 50)
    print(f"\nThe final calculated monoisotopic molecular weight is: {final_peptide_weight:.4f} Da")


if __name__ == '__main__':
    calculate_peptide_molecular_weight()