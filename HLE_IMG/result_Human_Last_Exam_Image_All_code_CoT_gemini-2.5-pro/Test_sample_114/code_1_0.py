import sys
import io

# Keep the original stdout
original_stdout = sys.stdout
# Redirect stdout to a string
sys.stdout = captured_output = io.StringIO()

def calculate_dou(formula_dict, formula_str):
    """Calculates and explains the Degree of Unsaturation."""
    c = formula_dict.get('C', 0)
    h = formula_dict.get('H', 0)
    n = formula_dict.get('N', 0)
    # The formula is: DoU = C - H/2 + N/2 + 1
    dou = c - (h / 2) + (n / 2) + 1
    print(f"For formula {formula_str} (C={c}, H={h}, N={n}):")
    # Using integer arithmetic for presentation if possible
    print(f"DoU = C - H/2 + N/2 + 1 = {c} - {h/2} + {n/2} + 1 = {int(dou)}")
    return int(dou)

# Molecular formulas from the problem description
product_A_formula = {'C': 9, 'H': 13, 'N': 1, 'O': 2}
product_B_formula = {'C': 10, 'H': 13, 'N': 1, 'O': 2}

print("### Rationale for Observed Regioselectivity ###\n")

print("1. Reaction Type and Intermediates:")
print("All three reactions are examples of 1,3-dipolar cycloadditions. The N-acyl amino acid starting materials react with acetic anhydride (Ac2O) and a base (TEA) to form a mesoionic 1,3-dipole known as a munchone. This munchone intermediate then reacts with the alkyne (methyl propiolate). The initial cycloadduct spontaneously loses CO2 to form the final pyrrole product.")
print("-" * 50)

print("2. Analysis of Reactions 1 & 2 (Product A):")
dou_A = calculate_dou(product_A_formula, "C9H13NO2")
print(f"The DoU of {dou_A} for Product A is consistent with a highly unsaturated structure containing a pyrrole ring (3 DU) and a carbonyl group (1 DU).")
print("\nReaction 2, the isotope labeling study, is crucial. It shows that the two methyl groups (one from the acetyl group, one from the alanine alpha-position) become scrambled. This implies that the munchone intermediate formed from the acyclic amino acid is not static. It exists in a rapid equilibrium that renders the two methyl groups chemically equivalent before the cycloaddition occurs. This scrambling likely proceeds through a reversible ring-opening to an acyclic ketene intermediate, which allows for rearrangement.")
print("-" * 50)

print("3. Analysis of Reaction 3 (Product B):")
dou_B = calculate_dou(product_B_formula, "C10H13NO2")
print(f"The DoU of {dou_B} for Product B is consistent with a bicyclic structure, containing a fused pyrrole-pyrrolidine ring system (4 DU) and a carbonyl group (1 DU).")
print("\nIn this case, the starting material is N-acetyl-proline. The munchone intermediate formed is fundamentally different because the proline's five-membered ring makes the intermediate a rigid, bicyclic system.")
print("-" * 50)

print("4. Conclusion: The Origin of Regioselectivity in Reaction 3")
print("The high regioselectivity observed in Reaction 3 is a direct consequence of its well-defined, conformationally rigid, and asymmetric intermediate.")
print("The munchone from proline has a methyl group at its C2 position and a methylene group (part of the fused ring) at its C4 position. These two positions are electronically and sterically distinct. Crucially, due to the rigid bicyclic structure, they cannot be scrambled by the rearrangement mechanism that occurs in Reaction 2.")
print("\nBecause the reaction proceeds through a single, structurally defined, asymmetric munchone, its subsequent cycloaddition with methyl propiolate is governed by strong and specific electronic preferences (as described by Frontier Molecular Orbital theory). This leads to the highly regioselective formation of only one regioisomer, Product B.")

# Restore stdout and print the captured output
sys.stdout = original_stdout
# We need to print the <<<answer>>> line at the end of all other prints.
# So we print the captured output first.
output = captured_output.getvalue()
print(output)
print("<<<The high regioselectivity in Reaction 3 results from a conformationally rigid and asymmetric munchone intermediate derived from proline, which, unlike the intermediate from Reaction 2, cannot undergo rearrangement to scramble its substituents before the cycloaddition occurs.>>>")