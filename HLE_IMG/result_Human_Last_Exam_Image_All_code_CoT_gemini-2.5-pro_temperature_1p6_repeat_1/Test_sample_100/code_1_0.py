import sys

def get_formula_string(atoms):
    """Generates a molecular formula string from a dictionary of atom counts."""
    # Standard order: C, H, then alphabetical for the rest.
    order = ['C', 'H'] + sorted([k for k in atoms if k not in ['C', 'H']])
    
    formula = ""
    for elem in order:
        count = atoms.get(elem, 0)
        if count > 0:
            formula += elem
            if count > 1:
                formula += str(count)
    return formula

def main():
    """
    Calculates the molecular formula of the product from the given reaction scheme.
    """
    # Base core molecule: 2-azabicyclo[2.2.1]hept-5-en-3-one (C6H7NO)
    core = {'C': 6, 'H': 7, 'N': 1, 'O': 1, 'F': 0}
    
    # Substituents
    # p-Methoxybenzyl group (C8H9O) replaces one H atom.
    pmb_group = {'C': 8, 'H': 9, 'O': 1}
    # Trifluoromethyl group (CF3) replaces one H atom.
    cf3_group = {'C': 1, 'F': 3}
    
    # --- Step 1: Starting Material ---
    c0 = core['C'] + pmb_group['C'] + cf3_group['C']
    h0 = core['H'] - 2 + pmb_group['H'] # -1H for PMB substitution, -1H for CF3 substitution
    n0 = core['N']
    o0 = core['O'] + pmb_group['O']
    f0 = core['F'] + cf3_group['F']
    start_formula = {'C': c0, 'H': h0, 'N': n0, 'O': o0, 'F': f0}

    print("Step 1: Determine the molecular formula of the starting material.")
    print(f"The core is C6H7NO. It is substituted with a PMB group (C8H9O) and a CF3 group (CF3).")
    print(f"The formula for the starting material is {get_formula_string(start_formula)}.")
    
    # --- Step 2: Intermediate 1 (CAN Deprotection) ---
    c1 = c0 - pmb_group['C']
    h1 = h0 - pmb_group['H'] + 1
    o1 = o0 - pmb_group['O']
    intermediate1_formula = {'C': c1, 'H': h1, 'N': n0, 'O': o1, 'F': f0}
    
    print("\nStep 2: Determine the molecular formula of Intermediate 1.")
    print("Reaction: Oxidative removal of the PMB group (loss of C8H9O, gain of H).")
    print(f"The formula for Intermediate 1 is {get_formula_string(intermediate1_formula)}.")

    # --- Step 3: Intermediate 2 (Hydrogenation) ---
    h2 = h1 + 2
    intermediate2_formula = {'C': c1, 'H': h2, 'N': n0, 'O': o1, 'F': f0}
    
    print("\nStep 3: Determine the molecular formula of Intermediate 2.")
    print("Reaction: Catalytic hydrogenation of the C=C double bond (addition of H2).")
    print(f"The formula for Intermediate 2 is {get_formula_string(intermediate2_formula)}.")

    # --- Step 4: Final Product (Hydrolysis) ---
    h3 = h2 + 2
    o3 = o1 + 1
    final_product_formula = {'C': c1, 'H': h3, 'N': n0, 'O': o3, 'F': f0}
    
    print("\nStep 4: Determine the molecular formula of the final product.")
    print("Reaction: Acid hydrolysis of the lactam (addition of H2O).")
    print(f"Final Formula Calculation:")
    print(f"C = {intermediate2_formula.get('C', 0)} = {final_product_formula.get('C', 0)}")
    print(f"H = {intermediate2_formula.get('H', 0)} + 2 = {final_product_formula.get('H', 0)}")
    print(f"F = {intermediate2_formula.get('F', 0)} = {final_product_formula.get('F', 0)}")
    print(f"N = {intermediate2_formula.get('N', 0)} = {final_product_formula.get('N', 0)}")
    print(f"O = {intermediate2_formula.get('O', 0)} + 1 = {final_product_formula.get('O', 0)}")

    final_formula_str = get_formula_string(final_product_formula)
    print(f"\nThe final molecular formula of the product is {final_formula_str}.")

if __name__ == "__main__":
    main()