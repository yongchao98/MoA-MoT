def solve_reaction():
    """
    Analyzes the given chemical reaction and provides details about the product.
    """

    # --- Reactant Information ---
    # The effective reactants are 2-aminopyridine, o-phthalaldehyde, and hydrogen cyanide (HCN),
    # as TMSCN reacts with water to form HCN and TMSOH.
    reactants = {
        "2-Aminopyridine": "C5H6N2",
        "o-Phthalaldehyde": "C8H6O2",
        "Hydrogen Cyanide": "HCN"
    }

    # --- Product Information ---
    # The reaction is a three-component synthesis leading to a specific heterocyclic product.
    # The overall reaction is: C5H6N2 + C8H6O2 + HCN -> C14H11N3O + H2O
    product_A = {
        "name": "3-cyano-2-(pyridin-2-yl)isoindolin-1-ol",
        "formula": "C14H11N3O"
    }
    
    water = {
        "name": "Water",
        "formula": "H2O"
    }

    # --- Atomic Weights (g/mol) ---
    atomic_weights = {
        'C': 12.011,
        'H': 1.008,
        'N': 14.007,
        'O': 15.999
    }

    def get_molecular_weight(formula):
        """Calculates the molecular weight of a compound from its formula."""
        import re
        
        # Find all atom-count pairs in the formula (e.g., C14, H11, N3, O)
        # re.findall will find patterns like 'C14', 'H11', 'N3', 'O1' (implicitly)
        elements = re.findall(r'([A-Z][a-z]*)(\d*)', formula)
        
        total_weight = 0.0
        for element, count in elements:
            # If count is an empty string, it means 1 atom
            num_atoms = int(count) if count else 1
            if element in atomic_weights:
                total_weight += num_atoms * atomic_weights[element]
        return total_weight

    # --- Output the Results ---
    print("--- Reaction Analysis ---")
    print("The reaction is a three-component synthesis of a substituted isoindoline.")
    print("\nReaction equation with stoichiometry and molecular formulas:")
    # The prompt requires outputting each number in the final equation.
    # Here, we print the balanced chemical equation.
    print("1 C5H6N2 + 1 C8H6O2 + 1 HCN -> 1 C14H11N3O + 1 H2O")

    print("\n--- Product A Details ---")
    print(f"Name: {product_A['name']}")
    print(f"Molecular Formula: {product_A['formula']}")
    
    # Calculate and print molecular weight
    mw = get_molecular_weight(product_A['formula'])
    print(f"Molecular Weight: {mw:.3f} g/mol")

solve_reaction()