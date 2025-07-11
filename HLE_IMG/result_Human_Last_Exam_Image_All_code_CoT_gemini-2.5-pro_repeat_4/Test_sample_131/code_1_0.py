def identify_alkaloid():
    """
    Identifies the alkaloid compound, provides its molecular formula,
    and calculates its molecular weight.
    """
    # Name of the compound
    compound_name = "Matrine"

    # Molecular formula: C15H24N2O
    atom_counts = {'C': 15, 'H': 24, 'N': 2, 'O': 1}

    # Standard atomic weights (g/mol)
    atomic_weights = {'C': 12.011, 'H': 1.008, 'N': 14.007, 'O': 15.999}

    # Print the name and formula
    print(f"The name of the alkaloid compound is: {compound_name}")
    print(f"Its molecular formula is C{atom_counts['C']}H{atom_counts['H']}N{atom_counts['N']}O{atom_counts['O']}.")
    print("\nCalculating the molecular weight:")

    # Calculate and build the equation string
    total_mw = 0
    equation_parts = []
    for atom, count in atom_counts.items():
        weight = atomic_weights[atom]
        total_mw += count * weight
        equation_parts.append(f"({count} * {weight})")

    equation_str = " + ".join(equation_parts)

    # Print the full equation and the result
    print(f"{equation_str} = {total_mw:.3f} g/mol")

if __name__ == "__main__":
    identify_alkaloid()