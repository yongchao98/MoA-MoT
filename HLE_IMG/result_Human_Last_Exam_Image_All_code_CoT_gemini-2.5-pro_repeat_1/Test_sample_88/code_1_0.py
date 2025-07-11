def solve_reaction_formulas():
    """
    Calculates and prints the derivation of product formulas based on
    proposed reaction pathways.
    """

    # Helper function to combine dictionaries representing molecular formulas
    def combine_formulas(f1, f2, sign=1):
        """Adds (sign=1) or subtracts (sign=-1) two formula dictionaries."""
        result = f1.copy()
        for atom, count in f2.items():
            result[atom] = result.get(atom, 0) + sign * count
        return result

    # Helper function to format and print a formula dictionary for display
    def print_formula_line(label, formula_dict):
        """Prints a labeled line with a formatted formula, aligning the columns."""
        atoms = ['C', 'H', 'N', 'O']
        # Format as "Atom count", e.g., "C 9  H 16 N 2  O 2"
        formula_str = " ".join([f"{atom} {formula_dict.get(atom, 0):<2}" for atom in atoms])
        print(f"{label:<47}: {formula_str}")

    print("This script verifies the molecular formulas of products A, B, and C.")
    print("The starting material (SM) is assumed to be 2-(1-pyrrolin-2-yl)proline, based on the drawing.")
    print("Its formula is C9 H16 N2 O2.\n")

    # --- Define base molecules and operations as dictionaries of atom counts ---
    SM = {'C': 9, 'H': 16, 'N': 2, 'O': 2}
    MP = {'C': 4, 'H': 4, 'O': 2} # Methyl Propiolate
    
    # Define operations as net changes in atom counts
    ACETYLATION = {'C': 2, 'H': 2, 'O': 1}  # Net change for N-H -> N-COCH3
    DECARBOXYLATION = {'C': 1, 'O': 2}     # Loss of CO2
    DEHYDROGENATION = {'H': 2}            # Loss of H2
    TRIPLE_DEHYDROGENATION = {'H': 6}     # Loss of 3*H2
    OXIDATION = {'O': 1}                  # Gain of one O atom

    # --- Calculation for Product C (C11H16N2O3) ---
    print("--- Derivation of Product C (Formula: C11 H16 N2 O3) ---")
    print_formula_line("Starting Material (SM)", SM)
    n_acetyl_sm = combine_formulas(SM, ACETYLATION)
    print_formula_line("+ Net effect of N-acetylation (+C2H2O)", ACETYLATION)
    print("-" * 70)
    print_formula_line("= N-acetyl-SM Intermediate", n_acetyl_sm)
    product_c = combine_formulas(n_acetyl_sm, DEHYDROGENATION, sign=-1)
    print_formula_line("- Dehydrogenation (loss of H2)", DEHYDROGENATION)
    print("-" * 70)
    print_formula_line("= Final Product C", product_c)
    print("Result matches the given formula for C.\n")

    # --- Calculation for Product A (C14H20N2O3) ---
    print("--- Derivation of Product A (Formula: C14 H20 N2 O3) ---")
    print_formula_line("N-acetyl-SM Intermediate (from above)", n_acetyl_sm)
    ylide_a = combine_formulas(n_acetyl_sm, DECARBOXYLATION, sign=-1)
    print_formula_line("- Decarboxylation (loss of CO2)", DECARBOXYLATION)
    print("-" * 70)
    print_formula_line("= Azomethine Ylide Intermediate", ylide_a)
    adduct_a = combine_formulas(ylide_a, MP)
    print_formula_line("+ Methyl Propiolate (C4H4O2)", MP)
    print("-" * 70)
    print_formula_line("= Initial Cycloadduct", adduct_a)
    product_a = combine_formulas(adduct_a, DEHYDROGENATION, sign=-1)
    print_formula_line("- Dehydrogenation (loss of H2)", DEHYDROGENATION)
    print("-" * 70)
    print_formula_line("= Final Product A", product_a)
    print("Result matches the given formula for A.\n")

    # --- Calculation for Product B (C12H14N2O3) ---
    print("--- Derivation of Product B (Formula: C12 H14 N2 O3) ---")
    print_formula_line("Starting Material (SM)", SM)
    ylide_b = combine_formulas(SM, DECARBOXYLATION, sign=-1)
    print_formula_line("- Decarboxylation (loss of CO2)", DECARBOXYLATION)
    print("-" * 70)
    print_formula_line("= Azomethine Ylide Intermediate", ylide_b)
    adduct_b = combine_formulas(ylide_b, MP)
    print_formula_line("+ Methyl Propiolate (C4H4O2)", MP)
    print("-" * 70)
    print_formula_line("= Initial Cycloadduct", adduct_b)
    oxidized_adduct_b = combine_formulas(adduct_b, TRIPLE_DEHYDROGENATION, sign=-1)
    print_formula_line("- 3x Dehydrogenation (loss of 3*H2)", TRIPLE_DEHYDROGENATION)
    product_b = combine_formulas(oxidized_adduct_b, OXIDATION)
    print_formula_line("+ Oxidation (gain of O)", OXIDATION)
    print("-" * 70)
    print_formula_line("= Final Product B", product_b)
    print("Result matches the given formula for B.\n")

# Run the calculation and print the results
solve_reaction_formulas()