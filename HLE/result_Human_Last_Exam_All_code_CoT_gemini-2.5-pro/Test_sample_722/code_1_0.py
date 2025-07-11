import sys

# Redirect print to a string to control output
from io import StringIO
old_stdout = sys.stdout
sys.stdout = mystdout = StringIO()

# --- Define Protein Masses ---
pA_theo = 25
pB_theo = 150
pC_theo = 60
pD_theo = 100
kinase = 40

print("### Analysis of SEC-MALS Experiments ###\n")

# --- Step 1: Analysis of Experiment 1 (Individual Proteins) ---
print("--- Step 1: Analysis of Experiment 1 (Individual Proteins) ---")
pA_exp1 = 50
pB_exp1 = 300
pC_exp1 = 60
pD_exp1 = 100

print(f"Protein A: Theoretical mass = {pA_theo} kDa, Experimental mass = {pA_exp1} kDa.")
oligomer_A = pA_exp1 / pA_theo
print(f"  Calculation: {pA_exp1} / {pA_theo} = {int(oligomer_A)}. Conclusion: Protein A is a homodimer.")
pA_native = pA_exp1

print(f"Protein B: Theoretical mass = {pB_theo} kDa, Experimental mass = {pB_exp1} kDa.")
oligomer_B = pB_exp1 / pB_theo
print(f"  Calculation: {pB_exp1} / {pB_theo} = {int(oligomer_B)}. Conclusion: Protein B is a homodimer.")
pB_native = pB_exp1

print(f"Protein C: Theoretical mass = {pC_theo} kDa, Experimental mass = {pC_exp1} kDa.")
oligomer_C = pC_exp1 / pC_theo
print(f"  Calculation: {pC_exp1} / {pC_theo} = {int(oligomer_C)}. Conclusion: Protein C is a monomer.")
pC_native = pC_exp1

print(f"Protein D: Theoretical mass = {pD_theo} kDa, Experimental mass = {pD_exp1} kDa.")
oligomer_D = pD_exp1 / pD_theo
print(f"  Calculation: {pD_exp1} / {pD_theo} = {int(oligomer_D)}. Conclusion: Protein D is a monomer.\n")
pD_native = pD_exp1

# --- Step 2: Analysis of Experiment 2 (Protein Mixture) ---
print("--- Step 2: Analysis of Experiment 2 (Protein Mixture) ---")
peak1_exp2 = 300
peak2_exp2 = 210
print(f"Observed peaks at {peak1_exp2} kDa and {peak2_exp2} kDa.")
print(f"The {peak1_exp2} kDa peak corresponds to the free Protein B dimer ({pB_native} kDa).")
acd_complex_mass = pA_native + pC_native + pD_native
print(f"The {peak2_exp2} kDa peak corresponds to a complex of Protein A dimer + Protein C + Protein D.")
print(f"  Equation: {pA_native} kDa (A-dimer) + {pC_native} kDa (C) + {pD_native} kDa (D) = {acd_complex_mass} kDa.")
print("Conclusion: In baseline conditions, Protein A-dimer, C, and D form a complex, while Protein B-dimer is free.\n")

# --- Step 3: Analysis of Experiment 3 (Mixture with Kinase) ---
print("--- Step 3: Analysis of Experiment 3 (Mixture with Kinase) ---")
peak1_exp3 = 25
peak2_exp3 = 40
peak3_exp3 = 460
print(f"Observed peaks at {peak1_exp3} kDa, {peak2_exp3} kDa, and {peak3_exp3} kDa.")
print(f"The {peak1_exp3} kDa peak is monomeric Protein A ({pA_theo} kDa). The dimer has dissociated due to phosphorylation.")
print(f"The {peak2_exp3} kDa peak is the free kinase ({kinase} kDa).")
bcd_complex_mass = pB_native + pC_native + pD_native
print(f"The {peak3_exp3} kDa peak is a complex of the remaining proteins: Protein B dimer + Protein C + Protein D.")
print(f"  Equation: {pB_native} kDa (B-dimer) + {pC_native} kDa (C) + {pD_native} kDa (D) = {bcd_complex_mass} kDa.")
print("Conclusion: Phosphorylation of Protein A decreases its affinity for C and D, causing it to dissociate. This allows B, C, and D to form a complex.\n")

# --- Step 4: Analysis of Experiment 4 (Dephosphorylation) ---
print("--- Step 4: Analysis of Experiment 4 (Dephosphorylation) ---")
peak1_exp4 = 50
peak2_exp4 = 460
print(f"Observed peaks at {peak1_exp4} kDa and {peak2_exp4} kDa.")
print(f"The {peak1_exp4} kDa peak is the re-formed Protein A dimer ({pA_native} kDa).")
print(f"The {peak2_exp4} kDa peak is the stable Protein B-C-D complex ({bcd_complex_mass} kDa).")
print("Conclusion: The B-C-D complex, once formed, is very stable and cannot be displaced by dephosphorylated Protein A dimer.\n")

# --- Step 5: Evaluation of Answer Choices ---
print("--- Step 5: Final Conclusion ---")
print("Based on the analysis:")
print(" - A, C, H are incorrect because there is no evidence Protein B is phosphorylated.")
print(" - B, D are incorrect because phosphorylation of Protein A *decreases* its affinity, it does not increase it.")
print(" - E is incorrect because Protein A is a monomer when phosphorylated (Exp 3), not *always* a dimer.")
print(" - F is incorrect because native Protein A is a dimer (Exp 1, 2, 4), not a monomer.")
print(" - G is incorrect because Exp 4 shows the B-C-D complex is stable against nonphosphorylated A, contradicting the affinity claim.")
print(" - I is incorrect because C and H are incorrect.")
print("\nAll provided options (A-I) contain statements that are contradicted by the experimental data.")

sys.stdout = old_stdout
print(mystdout.getvalue())
print("<<<J>>>")