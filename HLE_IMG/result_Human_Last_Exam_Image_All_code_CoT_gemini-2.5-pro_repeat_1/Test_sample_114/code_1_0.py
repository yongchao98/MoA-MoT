def calculate_dbe(formula):
    """
    Calculates the Degree of Unsaturation (Double Bond Equivalent)
    for a chemical formula string like 'C9H13NO2'.
    Handles C, H, O, N.
    """
    import re
    
    # Regex to find elements and their counts
    pattern = re.compile(r'([A-Z][a-z]?)(\d*)')
    
    # Default counts
    counts = {'C': 0, 'H': 0, 'N': 0}
    
    # Parse the formula string
    matches = pattern.findall(formula)
    for element, count in matches:
        if element in counts:
            if count == '':
                counts[element] += 1
            else:
                counts[element] += int(count)
    
    C = counts['C']
    H = counts['H']
    N = counts['N']
    
    # DBE formula: C - H/2 + N/2 + 1
    dbe = C - (H / 2) + (N / 2) + 1
    
    return int(dbe) if dbe == int(dbe) else dbe

# --- Calculations ---
formula_A = 'C9H13NO2'
formula_B = 'C10H13NO2'
dbe_A = calculate_dbe(formula_A)
dbe_B = calculate_dbe(formula_B)

# --- Output the Explanation ---
print("### Analysis of Product Structures ###")
print(f"The molecular formula for Product A is {formula_A}.")
print(f"The degree of unsaturation (DBE) for Product A is calculated as: C - H/2 + N/2 + 1 = 9 - 13/2 + 1/2 + 1 = {dbe_A}")
print(f"A DBE of {dbe_A} is consistent with a pyrrole ring structure (1 ring + 2 double bonds = 3 DBE) and an ester carbonyl group (1 DBE).")
print("-" * 40)
print(f"The molecular formula for Product B is {formula_B}.")
print(f"The degree of unsaturation (DBE) for Product B is calculated as: C - H/2 + N/2 + 1 = 10 - 13/2 + 1/2 + 1 = {dbe_B}")
print(f"A DBE of {dbe_B} is consistent with a bicyclic product containing two double bonds and one ester group (2 rings + 2 double bonds + 1 carbonyl = 5 DBE).")

print("\n### Rationale for Regioselectivity ###")
print("The key to understanding the regioselectivity lies in comparing the structures of the reactive intermediates, which are münchnones (mesoionic 1,3-dipoles), formed in each reaction.")
print("\n1. Lack of Selectivity in Reaction 2:")
print("   - The münchnone formed from N-acetyl-N-methyl-alanine in Reactions 1 and 2 has methyl groups at both of its reactive termini (C2 and C4).")
print("   - The 13C labeling experiment in Reaction 2 confirms that these two ends are so sterically and electronically similar that the incoming alkyne adds to either end with equal probability.")
print("   - This lack of preference results in no regioselectivity and a 1:1 mixture of products.")

print("\n2. High Selectivity in Reaction 3:")
print("   - The münchnone formed from N-acetyl-proline in Reaction 3 is a bicyclic system. Its C2 terminus has a methyl group, but its C4 terminus is a highly substituted bridgehead carbon, part of the fused pyrrolidine ring.")
print("   - This makes the two reactive termini of the dipole dramatically different, especially in terms of steric hindrance. The C4 position is far more crowded than the C2 position.")

print("\n3. Conclusion:")
print("   - The high regioselectivity observed in Reaction 3 is a direct result of this large steric difference. The cycloaddition proceeds via the transition state that minimizes steric repulsion.")
print("   - Consequently, the bulkier end of the alkyne (the -COOCH3 group) selectively attacks the less-hindered C2 terminus of the münchnone, leading to the formation of Product B as a single regioisomer.")
print("   - The results from Reactions 1 and 2 are critical because they provide the control experiment, demonstrating that without this significant, built-in steric bias, the reaction is unselective.")
