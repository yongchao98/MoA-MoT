def calculate_molecular_formula():
    """
    Calculates the molecular formula of the starting material from its IUPAC name.
    """
    # Define molecular formulas of the constituent parts
    # bicyclo[2.2.1]hept-5-ene: C7H10
    # At C7, CH2 is replaced by C(OMe)2. Net change from C7H10: -2H, +2C, +6H, +2O -> C9H14O2
    norbornene_ketal_core = {'C': 9, 'H': 14, 'O': 2}
    
    # At C2, 2H are replaced by -OH and the cyclopentenyl substituent
    # First, let's account for the -OH group. H is replaced by OH. Net change: +1O
    # Second, another H is replaced by the substituent.
    # So from the core, we subtract 2 H atoms and add the atoms of the two groups.
    base_skeleton = {
        'C': norbornene_ketal_core['C'],
        'H': norbornene_ketal_core['H'] - 2,
        'O': norbornene_ketal_core['O']
    }
    
    oh_group = {'H': 1, 'O': 1}
    
    # Substituent: (S)-4-((tert-butyldimethylsilyl)oxy)cyclopent-1-en-1-yl
    # cyclopent-1-en-1-yl group: C5H7
    # (tert-butyldimethylsilyl)oxy group (OTBS): O-Si(Me)2(tBu) -> OSiC6H15
    # The substituent is C5H6(OTBS). Formula = C(5+6) H(6+15) O(1) Si(1) -> C11H21OSi
    substituent = {'C': 11, 'H': 21, 'O': 1, 'Si': 1}

    # Sum all parts for the final formula
    final_formula = {}
    all_components = [base_skeleton, oh_group, substituent]
    for comp in all_components:
        for element, count in comp.items():
            final_formula[element] = final_formula.get(element, 0) + count
            
    return final_formula

def format_formula(formula_dict):
    """Formats a dictionary into a chemical formula string."""
    # Order: C, H, then alphabetical for others
    order = ['C', 'H', 'O', 'Si']
    parts = []
    for element in order:
        if element in formula_dict:
            parts.append(f"{element}{formula_dict[element]}")
    return " ".join(parts)

# Main execution
reactant_name = "(1S,2R,4S)-2-((S)-4-((tert-butyldimethylsilyl)oxy)cyclopent-1-en-1-yl)-7,7-dimethoxybicyclo[2.2.1]hept-5-en-2-ol"
product_description = "A fused bicyclic α,β-unsaturated ketone with a 9-membered and a 5-membered ring."

# The reaction is an Anionic Oxy-Cope rearrangement, which preserves the molecular formula.
formula_dict = calculate_molecular_formula()
formula_str = format_formula(formula_dict)

print("Reaction Analysis:")
print("------------------")
print("The reaction is an Anionic Oxy-Cope [3,3]-sigmatropic rearrangement.")
print("\nFinal Equation:")
print("--------------")
print(f"Reactant: {reactant_name}")
print(f"Product: {product_description}")
print("\nMolecular Formula of Reactant and Product:")
print("-----------------------------------------")
# "output each number in the final equation!"
# Printing the counts for each element in the final formula.
for element, count in sorted(formula_dict.items()):
    print(f"Number of {element} atoms: {count}")
print(f"\nFormula: {formula_str}")

<<<A fused bicyclic α,β-unsaturated ketone with a 9-membered and a 5-membered ring.>>>