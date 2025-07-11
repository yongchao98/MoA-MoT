def calculate_dbe(C, H, N, X=0):
    """Calculates the degree of unsaturation for a molecule with formula CxHyNzXx."""
    return C - (H / 2) + (N / 2) + 1

# --- Data from the problem ---
# Product A
formula_A = {'C': 9, 'H': 13, 'N': 1, 'O': 2}
dbe_A = calculate_dbe(formula_A['C'], formula_A['H'], formula_A['N'])

# Product B
formula_B = {'C': 10, 'H': 13, 'N': 1, 'O': 2}
dbe_B = calculate_dbe(formula_B['C'], formula_B['H'], formula_B['N'])

# --- Analysis and Explanation ---

print("--- Analysis of Products and Reaction Mechanism ---")

# Step 1: Analyze Product A
print("\n[Product A Analysis]")
print(f"Molecular Formula (from HRMS): C{formula_A['C']}H{formula_A['H']}N{formula_A['N']}O{formula_A['O']}")
print(f"Degree of Unsaturation (DBE): {dbe_A}")
print("A DBE of 4.0 is consistent with a pyrrole ring (1 ring + 2 double bonds = 3 DBEs) and a carbonyl group (1 DBE).")
print("This, along with the NMR data, confirms Product A is a substituted pyrrole.")

# Step 2: Analyze Product B
print("\n[Product B Analysis]")
print(f"Molecular Formula (from HRMS): C{formula_B['C']}H{formula_B['H']}N{formula_B['N']}O{formula_B['O']}")
print(f"Degree of Unsaturation (DBE): {dbe_B}")
print("A DBE of 5.0 is consistent with a fused bicyclic structure containing a pyrrole ring and a carbonyl group (e.g., a dihydroindolizine derivative with 2 rings + 2 double bonds = 4 DBEs, plus a carbonyl = 5 DBEs).")

# Step 3: Explain the Regioselectivity
print("\n--- Explanation of Regioselectivity ---")
print("\nThe reaction is a 1,3-dipolar cycloaddition. The N-acyl amino acid reacts with Ac2O and cat. TEA to form a mesoionic intermediate (a 'munchone'), which then reacts with the alkyne.")

print("\n[Reactions 1 and 2: Acyclic Case & The Role of Labeling]")
print("In Reaction 1, the acyclic amino acid reacts under the specified conditions (cat. TEA, Ac2O, 65°C, 3 hr) to form an acyclic munchone intermediate.")
print("Reaction 2 is the key to understanding the system. The observation that 13C labels get scrambled to give a 1:1 ratio of products means the two precursor methyl groups (from the acetyl and alanine parts) become equivalent before the final product is formed.")
print("This scrambling occurs because the flexible, acyclic munchone intermediate undergoes a rapid [3,3]-sigmatropic rearrangement. This process interconverts two isomeric forms of the munchone, scrambling the positions of its substituents. This dynamic equilibrium is the reason for the lack of regioselectivity.")

print("\n[Reaction 3: Cyclic Case & The Origin of Selectivity]")
print("In Reaction 3, the cyclic amino acid (N-acetylproline) forms a rigid, bicyclic munchone intermediate.")
print("Crucially, this intermediate's fused-ring structure makes it unable to undergo the [3,3]-sigmatropic rearrangement. The atoms are sterically 'locked' in place.")
print("Because the intermediate is configurationally stable and cannot equilibrate, its reaction with the 5 equivalents of methyl propiolate is highly regioselective.")
print("This results in the formation of a single regioisomeric product, Product B.")

print("\n[Conclusion]")
print("The high regioselectivity in Reaction 3 is caused by its rigid, bicyclic munchone intermediate, which prevents the substituent-scrambling rearrangement seen in the acyclic case. The results from Reactions 1 and 2 are vital as they demonstrate this scrambling phenomenon in a flexible system, providing the necessary contrast to explain the selectivity of the rigid system.")