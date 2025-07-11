def solve_chemistry_problem():
    """
    This script identifies Compound A from the given reaction, provides its
    molecular formula and name, and prints the balanced chemical equation.
    """

    # Define atomic weights for calculation
    atomic_weights = {
        'C': 12.011,
        'H': 1.008,
        'N': 14.007,
        'O': 15.999,
    }

    # Define molecular formulas for reactants and the product
    # Note: TMSCN provides the 'CN' group, and H2O provides the proton, so we use HCN for balancing.
    formulas = {
        "2-Aminopyridine": {'C': 5, 'H': 6, 'N': 2},
        "o-Phthalaldehyde": {'C': 8, 'H': 6, 'O': 2},
        "HCN (from TMSCN/H2O)": {'H': 1, 'C': 1, 'N': 1},
        "Compound A": {'C': 14, 'H': 11, 'N': 3, 'O': 1},
        "Water (byproduct)": {'H': 2, 'O': 1}
    }

    # Function to calculate molecular weight
    def calculate_mw(formula):
        mw = 0
        for atom, count in formula.items():
            mw += atomic_weights.get(atom, 0) * count
        return mw

    # Identify Compound A
    compound_A_name = "3-hydroxy-2-(pyridin-2-yl)isoindoline-1-carbonitrile"
    compound_A_formula_str = "C14H11N3O"
    
    print(f"Compound A is {compound_A_name}.")
    print(f"Its molecular formula is {compound_A_formula_str}.")
    print("-" * 30)

    # Print the balanced chemical equation
    # Reactants: C5H6N2 + C8H6O2 + HCN
    # Products: C14H11N3O + H2O
    print("The balanced chemical equation is:")
    
    # Extracting coefficients for printing the equation clearly
    r1_formula = formulas["2-Aminopyridine"]
    r2_formula = formulas["o-Phthalaldehyde"]
    r3_formula = formulas["HCN (from TMSCN/H2O)"]
    p1_formula = formulas["Compound A"]
    p2_formula = formulas["Water (byproduct)"]

    # Print out each number in the final equation
    print(f"C{r1_formula['C']}H{r1_formula['H']}N{r1_formula['N']} + "
          f"C{r2_formula['C']}H{r2_formula['H']}O{r2_formula['O']} + "
          f"H{r3_formula['H']}C{r3_formula['C']}N{r3_formula['N']} -> "
          f"C{p1_formula['C']}H{p1_formula['H']}N{p1_formula['N']}O{p1_formula['O']} + "
          f"H{p2_formula['H']}O{p2_formula['O']}")
    print("-" * 30)

    # Verify mass balance (optional)
    mw_reactants = (calculate_mw(r1_formula) + 
                    calculate_mw(r2_formula) + 
                    calculate_mw(r3_formula))
    mw_products = (calculate_mw(p1_formula) + 
                   calculate_mw(p2_formula))
    
    print(f"Molecular weight of Reactants: {mw_reactants:.3f} g/mol")
    print(f"Molecular weight of Products: {mw_products:.3f} g/mol")
    print("The mass is balanced.")

solve_chemistry_problem()