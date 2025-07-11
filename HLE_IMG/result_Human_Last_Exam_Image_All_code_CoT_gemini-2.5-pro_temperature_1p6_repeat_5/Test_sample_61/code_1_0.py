#
# This script calculates the molecular formula of compound A based on the reaction provided.
# The chemical logic is embedded in the calculation steps.
#

def calculate_molecular_formula():
    """
    Calculates the molecular formula of the product A by simulating the reaction
    using atomic arithmetic.
    """
    # 1. Define the atomic composition of the reactants and key functional groups.
    # Compound 1 is Methyl 1-oxo-2,3-dihydro-1H-indene-2-carboxylate.
    # C: 6 (benzene ring) + 3 (cyclopentane ring) + 1 (ester C=O) + 1 (ester O-Me) = 11
    # H: 4 (aromatic) + 1 (alpha-H) + 2 (beta-CH2) + 3 (methyl) = 10
    # O: 1 (ketone) + 2 (ester) = 3
    compound1 = {'C': 11, 'H': 10, 'O': 3}

    # The benzyl group from benzyl bromide (C6H5CH2-).
    benzyl_group = {'C': 7, 'H': 7}

    # The methoxycarbonyl group (-COOMe) that is removed during saponification and decarboxylation.
    coome_group = {'C': 2, 'H': 3, 'O': 2}

    print("Step-by-step calculation of the molecular formula for compound A:")
    print("-" * 60)

    # --- Step 1: Alkylation ---
    # In the alkylation step, one hydrogen atom on compound 1 is replaced by the benzyl group.
    intermediate_C = compound1['C'] - 0 + benzyl_group['C']
    intermediate_H = compound1['H'] - 1 + benzyl_group['H']
    intermediate_O = compound1['O'] - 0 + 0

    print("Step 1: Alkylation")
    print("An H atom on compound 1 is replaced by a benzyl group.")
    print(f"Formula of Compound 1: C{compound1['C']} H{compound1['H']} O{compound1['O']}")
    print(f"Formula of Benzyl group: C{benzyl_group['C']} H{benzyl_group['H']}")
    print(f"Equation for C atoms: {compound1['C']} + {benzyl_group['C']} = {intermediate_C}")
    print(f"Equation for H atoms: {compound1['H']} - 1 + {benzyl_group['H']} = {intermediate_H}")
    print(f"Equation for O atoms: {compound1['O']} + 0 = {intermediate_O}")
    print(f"--> Formula after alkylation: C{intermediate_C}H{intermediate_H}O{intermediate_O}\n")


    # --- Step 2: Saponification and Decarboxylation ---
    # This combined process results in the net replacement of the -COOMe group with an H atom.
    final_C = intermediate_C - coome_group['C']
    final_H = intermediate_H - coome_group['H'] + 1
    final_O = intermediate_O - coome_group['O']

    print("Step 2: Saponification & Decarboxylation")
    print("The -COOMe group is replaced by an H atom.")
    print(f"Formula of intermediate: C{intermediate_C} H{intermediate_H} O{intermediate_O}")
    print(f"Formula of -COOMe group: C{coome_group['C']} H{coome_group['H']} O{coome_group['O']}")
    print(f"Equation for C atoms: {intermediate_C} - {coome_group['C']} = {final_C}")
    print(f"Equation for H atoms: {intermediate_H} - {coome_group['H']} + 1 = {final_H}")
    print(f"Equation for O atoms: {intermediate_O} - {coome_group['O']} = {final_O}")

    # --- Final Result ---
    # Construct the final formula string, handling the case where O=1.
    final_formula_str = f"C{final_C}H{final_H}"
    if final_O == 1:
        final_formula_str += "O"
    elif final_O > 1:
        final_formula_str += f"O{final_O}"
    
    print("-" * 60)
    print(f"The final molecular formula of compound A is: {final_formula_str}")

# Run the calculation
calculate_molecular_formula()