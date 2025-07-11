def calculate_product_formula():
    """
    Calculates the molecular formula of the final product through a three-step synthesis.
    """
    
    # Step 0: Determine the molecular formula of the starting material.
    # The starting material is 2-(4-methoxybenzyl)-5-(trifluoromethyl)-2-azabicyclo[2.2.1]hept-5-en-3-one.
    # Bicyclic core (C6H5F3NO skeleton): 6 C, 5 H, 3 F, 1 N, 1 O
    # PMB group (p-methoxybenzyl, -CH2-C6H4-OCH3): 8 C, 9 H, 1 O
    # Total Starting Material: C(6+8) H(5+9) F3 N1 O(1+1) = C14H14F3NO2
    formula = {'C': 14, 'H': 14, 'F': 3, 'N': 1, 'O': 2}
    
    print("Step-by-step determination of the molecular formula:")
    print("--------------------------------------------------")
    print(f"Starting Material Formula: C{formula['C']}H{formula['H']}F{formula['F']}N{formula['N']}O{formula['O']}")
    
    # Step 1: Deprotection of the PMB group with CAN.
    # The PMB group (C8H9O) is removed and replaced with a single Hydrogen atom on the Nitrogen.
    # Net change: -8 C, -8 H, -1 O.
    formula['C'] -= 8
    formula['H'] -= 9
    formula['O'] -= 1
    formula['H'] += 1  # Add H to the nitrogen
    
    print(f"Reaction 1 (CAN, ACN/H2O): PMB group removal.")
    print(f"Intermediate 1 Formula: C{formula['C']}H{formula['H']}F{formula['F']}N{formula['N']}O{formula['O']}")
    
    # Step 2: Hydrogenation of the C=C double bond with Pd/C, H2.
    # This adds two Hydrogen atoms across the double bond.
    formula['H'] += 2
    
    print(f"Reaction 2 (Pd/C, H2): Hydrogenation of the alkene.")
    print(f"Intermediate 2 Formula: C{formula['C']}H{formula['H']}F{formula['F']}N{formula['N']}O{formula['O']}")

    # Step 3: Hydrolysis of the lactam (cyclic amide) with HCl.
    # This adds one molecule of water (H2O), breaking the C-N amide bond.
    formula['H'] += 2
    formula['O'] += 1
    
    print(f"Reaction 3 (4 N HCl, 70 C): Lactam hydrolysis.")
    final_formula_str = f"C{formula['C']}H{formula['H']}F{formula['F']}N{formula['N']}O{formula['O']}"
    
    # Format the final formula, omitting 1s for N and F3 instead of F3.
    final_formula_str_formatted = f"C{formula['C']}H{formula['H']}F{formula['F']}NO{formula['O']}"

    print(f"Final Product Formula: {final_formula_str_formatted}")

calculate_product_formula()
