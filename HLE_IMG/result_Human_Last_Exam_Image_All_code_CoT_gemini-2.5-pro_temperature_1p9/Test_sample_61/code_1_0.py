def solve_chemistry_problem():
    """
    This script calculates the molecular formula of compound A based on the given reaction.
    """
    # Step 1: Define the molecular formula of the reactants.
    # Compound 1: methyl 1-oxo-2,3-dihydro-1H-indene-2-carboxylate
    # Structure breakdown: C6H4 (benzene ring) + C=O (C1) + CH (C2) + CH2 (C3) + COOMe (ester)
    # C: 6 + 1 + 1 + 1 + 2 = 11
    # H: 4 + 1 + 2 + 3 = 10
    # O: 1 + 2 = 3
    comp1 = {'C': 11, 'H': 10, 'O': 3}

    # Benzyl group from Compound 2 (benzyl bromide, C7H7Br)
    # The reacting part is the benzyl group C7H7.
    benzyl_group = {'C': 7, 'H': 7}

    print("Step-by-step calculation for the molecular formula of Compound A:")

    print("\n1. Molecular formula of Reactant 1:")
    print(f"   C: {comp1['C']}, H: {comp1['H']}, O: {comp1['O']} => C{comp1['C']}H{comp1['H']}O{comp1['O']}")

    # Step 2: Calculate the formula of the alkylated intermediate.
    # The reaction replaces one hydrogen atom from Compound 1 with a benzyl group.
    inter_c = comp1['C'] + benzyl_group['C']
    inter_h = comp1['H'] - 1 + benzyl_group['H']
    inter_o = comp1['O']

    print("\n2. After alkylation (replacing H with a benzyl group C7H7):")
    print("   The formula becomes C(11+7)H(10-1+7)O3")
    print(f"   C: {comp1['C']} + {benzyl_group['C']} = {inter_c}")
    print(f"   H: {comp1['H']} - 1 + {benzyl_group['H']} = {inter_h}")
    print(f"   O: {comp1['O']} = {inter_o}")
    print(f"   Intermediate formula: C{inter_c}H{inter_h}O{inter_o}")


    # Step 3: Calculate the formula of the final product A after saponification.
    # The methyl ester group (-COOCH3) is hydrolyzed to a carboxylic acid (-COOH).
    # This corresponds to replacing a methyl group (-CH3) with a hydrogen atom (-H).
    # Net change is a loss of one C and two H atoms (-CH2).
    final_a_c = inter_c - 1
    final_a_h = inter_h - 2
    final_a_o = inter_o

    print("\n3. After saponification (hydrolysis of methyl ester, -CH3 replaced by -H):")
    print("   The formula changes by a net loss of CH2.")
    print("   The final formula is C(18-1)H(16-2)O3")
    print(f"   C: {inter_c} - 1 = {final_a_c}")
    print(f"   H: {inter_h} - 2 = {final_a_h}")
    print(f"   O: {inter_o} = {final_a_o}")


    # Step 4: Display the final result.
    final_formula = f"C{final_a_c}H{final_a_h}O{final_a_o}"
    print("\nTherefore, the final molecular formula of compound A is:")
    print(final_formula)

solve_chemistry_problem()