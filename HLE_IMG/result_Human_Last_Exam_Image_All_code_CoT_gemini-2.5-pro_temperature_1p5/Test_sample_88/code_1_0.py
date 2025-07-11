def print_formula(label, atoms):
    """Helper function to print a chemical formula from a dictionary of atoms."""
    return f"{label}(C{atoms.get('C', 0)}H{atoms.get('H', 0)}N{atoms.get('N', 0)}O{atoms.get('O', 0)})"

def main():
    """
    Calculates and verifies the molecular formulas of products A, B, and C
    based on hypothesized reaction pathways.
    """
    # Molecular formulas of reactants and known molecules
    sm = {'C': 9, 'H': 14, 'N': 2, 'O': 2}      # Starting Material (SM)
    mp = {'C': 4, 'H': 4, 'O': 2}        # Methyl Propiolate (MP)
    acetyl_group = {'C': 2, 'H': 3, 'O': 1} # Acetyl group from Acetic Anhydride
    H_atom = {'H': 1}                   # A single hydrogen atom
    co2 = {'C': 1, 'O': 2}              # Carbon Dioxide
    meoh = {'C': 1, 'H': 4, 'O': 1}     # Methanol

    # Known molecular formulas of the products
    product_A = {'C': 14, 'H': 20, 'N': 2, 'O': 3}
    product_B = {'C': 12, 'H': 14, 'N': 2, 'O': 3}
    product_C = {'C': 11, 'H': 16, 'N': 2, 'O': 3}

    print("--- Analysis of Product Formation ---")

    # --- Pathway for Product C ---
    print("\n1. Pathway for Product C: N-Acetylation")
    print("Product C is formed by the N-acetylation of the starting material (SM).")
    c_calc = {
        'C': sm['C'] + acetyl_group['C'],
        'H': sm['H'] + acetyl_group['H'] - H_atom['H'],
        'N': sm['N'],
        'O': sm['O'] + acetyl_group['O']
    }
    print(f"   Equation: SM + Acetyl_Group - H = Product_C")
    print(f"   Calculation: C({sm['C']}+{acetyl_group['C']}) H({sm['H']}+{acetyl_group['H']}-{H_atom['H']}) N({sm['N']}) O({sm['O']}+{acetyl_group['O']})")
    print(f"   Result: {print_formula('C_calc', c_calc)}")
    print(f"   Expected: {print_formula('Product C', product_C)}")
    if c_calc == product_C:
        print("   Verification: SUCCESS, calculated formula matches Product C.")
    else:
        print("   Verification: FAILED.")


    # --- Pathway for Product A ---
    print("\n2. Pathway for Product A: 1,3-Dipolar Cycloaddition")
    print("Product A is formed from Product C via an azomethine ylide intermediate.")
    
    # Step 2a: Formation of the Azomethine Ylide from Product C
    print("   Step 2a: Product C undergoes thermal decarboxylation to form an azomethine ylide.")
    ylide_calc = {
        'C': product_C['C'] - co2['C'],
        'H': product_C['H'],
        'N': product_C['N'],
        'O': product_C['O'] - co2['O']
    }
    print(f"      Equation: Product_C - CO2 = Ylide")
    print(f"      Calculation: C({product_C['C']}-{co2['C']}) H({product_C['H']}) N({product_C['N']}) O({product_C['O']}-{co2['O']})")
    print(f"      Result: {print_formula('Ylide', ylide_calc)}")

    # Step 2b: Cycloaddition of the Ylide with Methyl Propiolate (MP)
    print("\n   Step 2b: The ylide reacts with methyl propiolate (MP) in a [3+2] cycloaddition.")
    a_calc = {
        'C': ylide_calc['C'] + mp['C'],
        'H': ylide_calc['H'] + mp['H'],
        'N': ylide_calc['N'],
        'O': ylide_calc['O'] + mp['O']
    }
    print(f"      Equation: Ylide + MP = Product_A")
    print(f"      Calculation: C({ylide_calc['C']}+{mp['C']}) H({ylide_calc['H']}+{mp['H']}) N({ylide_calc['N']}) O({ylide_calc['O']}+{mp['O']})")
    print(f"      Result: {print_formula('A_calc', a_calc)}")
    print(f"      Expected: {print_formula('Product A', product_A)}")
    if a_calc == product_A:
        print("      Verification: SUCCESS, calculated formula matches Product A.")
    else:
        print("      Verification: FAILED.")

    # --- Pathway for Product B ---
    print("\n3. Pathway for Product B: Michael Addition and Condensation")
    print("Product B is formed by the reaction of SM with MP, followed by the elimination of methanol.")
    b_calc = {
        'C': sm['C'] + mp['C'] - meoh['C'],
        'H': sm['H'] + mp['H'] - meoh['H'],
        'N': sm['N'],
        'O': sm['O'] + mp['O'] - meoh['O']
    }
    print(f"   Equation: SM + MP - Methanol = Product_B")
    print(f"   Calculation: C({sm['C']}+{mp['C']}-{meoh['C']}) H({sm['H']}+{mp['H']}-{meoh['H']}) N({sm['N']}) O({sm['O']}+{mp['O']}-{meoh['O']})")
    print(f"   Result: {print_formula('B_calc', b_calc)}")
    print(f"   Expected: {print_formula('Product B', product_B)}")
    if b_calc == product_B:
        print("   Verification: SUCCESS, calculated formula matches Product B.")
    else:
        print("   Verification: FAILED.")

if __name__ == "__main__":
    main()