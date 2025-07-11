def format_formula(atom_counts):
    """Formats a dictionary of atom counts into a molecular formula string."""
    formula_str = ""
    # Standard order: C, H, then alphabetical for the rest
    order = ['C', 'H', 'F', 'N', 'O']
    for atom in order:
        if atom in atom_counts and atom_counts[atom] > 0:
            formula_str += atom
            count = atom_counts[atom]
            if count > 1:
                formula_str += str(count)
    return formula_str

def main():
    """
    Calculates the molecular formula of the product through the reaction series.
    """
    print("Determining the molecular formula of the final product step-by-step:\n")

    # Step 1: Starting Material
    # Structure: 2-azabicyclo[2.2.1]hept-5-en-3-one core with a CF3 group and a PMB group.
    # Core (C6H5F3NO) + PMB group (C8H9O) = C15H14F3NO2
    start_formula = {'C': 15, 'H': 14, 'F': 3, 'N': 1, 'O': 2}
    print(f"1. Starting Material Formula: {format_formula(start_formula)}")

    # Step 2: Reaction 1 (CAN deprotection) -> Intermediate 1
    # Removes the p-methoxybenzyl (PMB) group (C8H9O) and adds a hydrogen (H).
    intermediate1_formula = start_formula.copy()
    intermediate1_formula['C'] -= 8
    intermediate1_formula['H'] -= 9
    intermediate1_formula['O'] -= 1
    intermediate1_formula['H'] += 1
    print(f"2. Intermediate 1 Formula (after deprotection): {format_formula(intermediate1_formula)}")

    # Step 3: Reaction 2 (Hydrogenation) -> Intermediate 2
    # Reduces the C=C double bond by adding H2.
    intermediate2_formula = intermediate1_formula.copy()
    intermediate2_formula['H'] += 2
    print(f"3. Intermediate 2 Formula (after hydrogenation): {format_formula(intermediate2_formula)}")

    # Step 4: Reaction 3 (Hydrolysis) -> Product
    # Hydrolyzes the lactam (cyclic amide) by adding H2O.
    product_formula = intermediate2_formula.copy()
    product_formula['H'] += 2
    product_formula['O'] += 1
    product_str = format_formula(product_formula)
    print(f"4. Final Product Formula (after hydrolysis): {product_str}")
    
    print("\n--- Final Calculation Breakdown ---")
    inter2_str = format_formula(intermediate2_formula)
    h2o_str = format_formula({'H': 2, 'O': 1})
    
    print(f"The final reaction is the hydrolysis of Intermediate 2:")
    print(f"{inter2_str} + {h2o_str} -> {product_str}")

    print("\nThe molecular formula of the final product is composed of:")
    print(f"Carbon (C) atoms: {product_formula['C']}")
    print(f"Hydrogen (H) atoms: {product_formula['H']}")
    print(f"Fluorine (F) atoms: {product_formula['F']}")
    print(f"Nitrogen (N) atoms: {product_formula['N']}")
    print(f"Oxygen (O) atoms: {product_formula['O']}")

if __name__ == "__main__":
    main()