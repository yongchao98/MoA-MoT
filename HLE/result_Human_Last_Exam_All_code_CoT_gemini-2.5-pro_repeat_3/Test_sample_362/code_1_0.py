def get_product_of_wittig_reaction():
    """
    Calculates and displays the products of a Wittig reaction between
    pivalaldehyde and (2-(2-chlorophenyl)ethylidene)triphenylphosphane.
    """
    # Molecular formula is represented as a dictionary of {atom: count}

    # 1. Define Reactant Formulas
    # Pivalaldehyde: (CH3)3C-CHO -> C5H10O
    pivalaldehyde_formula = {'C': 5, 'H': 10, 'O': 1}

    # Ylide: PPh3=CH-CH2-(C6H4Cl)
    # PPh3 part: P(C6H5)3 -> C18 H15 P
    # Ylide carbon chain part: =CH-CH2-(C6H4Cl) -> C8 H7 Cl
    ylide_formula = {'C': 18 + 8, 'H': 15 + 7, 'Cl': 1, 'P': 1}

    # 2. Determine Product Formulas based on the Wittig reaction mechanism
    # The reaction swaps the aldehyde's =O with the ylide's =C< group.

    # Main Alkene Product = (pivalaldehyde - O) + (ylide - PPh3)
    aldehyde_fragment = {'C': 5, 'H': 10}
    ylide_fragment = {'C': 8, 'H': 7, 'Cl': 1}
    
    alkene_product_formula = {}
    all_atoms = set(aldehyde_fragment.keys()) | set(ylide_fragment.keys())
    for atom in sorted(list(all_atoms)):
        alkene_product_formula[atom] = aldehyde_fragment.get(atom, 0) + ylide_fragment.get(atom, 0)

    # Byproduct = PPh3 + O
    phosphine_oxide_formula = {'C': 18, 'H': 15, 'P': 1, 'O': 1}

    # Helper function to format the formula string as requested
    def format_formula(formula_dict):
        """Formats a formula dict into a string like C(x)H(y)..."""
        parts = []
        # Standard order: C, H, then alphabetical for the rest
        order = ['C', 'H'] + sorted([k for k in formula_dict if k not in ['C', 'H']])
        for atom in order:
            if atom in formula_dict:
                count = formula_dict[atom]
                # Avoid printing (1) for clarity, unless it's the only atom
                if count > 1 or len(formula_dict) == 1:
                    parts.append(f"{atom}({count})")
                else:
                    parts.append(atom)
        return "".join(parts)

    # 3. Print the results
    product_name = "(Z)-1-(2-chlorophenyl)-4,4-dimethylpent-2-ene"
    byproduct_name = "Triphenylphosphine oxide"

    print("The Wittig reaction combines an aldehyde and a ylide to form an alkene and a phosphine oxide.")
    print("\n--- Reactants ---")
    print(f"Aldehyde: Pivalaldehyde")
    print(f"Ylide: (2-(2-chlorophenyl)ethylidene)triphenylphosphane")
    
    print("\n--- Products ---")
    print(f"Major Alkene Product: {product_name}")
    print(f"Byproduct: {byproduct_name}")

    print("\n--- Balanced Chemical Equation ---")
    pivalaldehyde_str = format_formula(pivalaldehyde_formula)
    ylide_str = format_formula(ylide_formula)
    alkene_str = format_formula(alkene_product_formula)
    phosphine_oxide_str = format_formula(phosphine_oxide_formula)
    
    # Final equation with atom counts
    print(f"The final equation with atom counts is:")
    print(f"{pivalaldehyde_str} + {ylide_str} -> {alkene_str} + {phosphine_oxide_str}")


get_product_of_wittig_reaction()