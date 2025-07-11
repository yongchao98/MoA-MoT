import sys
import io

# Capture the original stdout to restore it later if needed
original_stdout = sys.stdout
# Create a string buffer to hold the output
output_buffer = io.StringIO()
# Redirect stdout to the buffer
sys.stdout = output_buffer

# --- Define Protein Masses and Experimental Data ---
protein_A_theo_mass = 25
protein_B_theo_mass = 150
protein_C_theo_mass = 60
protein_D_theo_mass = 100
kinase_mass = 40

# --- Step-by-Step Analysis ---

print("Step-by-step Analysis of the SEC-MALS Experiments\n")

# 1. Analyze Experiment 1
print("--- Analysis of Experiment 1 (Individual Proteins) ---")
exp1_A_mass = 50
exp1_B_mass = 300
exp1_C_mass = 60
exp1_D_mass = 100
A_oligomer_factor = exp1_A_mass / protein_A_theo_mass
B_oligomer_factor = exp1_B_mass / protein_B_theo_mass
A_dimer_mass = exp1_A_mass
B_dimer_mass = exp1_B_mass
C_monomer_mass = exp1_C_mass
D_monomer_mass = exp1_D_mass

print(f"Protein A: Runs at {exp1_A_mass} kDa (theoretical: {protein_A_theo_mass} kDa). This is a homodimer ({int(A_oligomer_factor)} * {protein_A_theo_mass} = {exp1_A_mass}).")
print(f"Protein B: Runs at {exp1_B_mass} kDa (theoretical: {protein_B_theo_mass} kDa). This is a homodimer ({int(B_oligomer_factor)} * {protein_B_theo_mass} = {exp1_B_mass}).")
print(f"Protein C: Runs at {exp1_C_mass} kDa (theoretical: {protein_C_theo_mass} kDa). This is a monomer.")
print(f"Protein D: Runs at {exp1_D_mass} kDa (theoretical: {protein_D_theo_mass} kDa). This is a monomer.")
print("Conclusion 1: Proteins A and B exist as stable homodimers, while C and D are monomers.\n")

# 2. Analyze Experiment 2
print("--- Analysis of Experiment 2 (Mixture without Kinase) ---")
exp2_peak1 = 300
exp2_peak2 = 210
complex_ACD_mass = A_dimer_mass + C_monomer_mass + D_monomer_mass
print(f"Peak 1 at {exp2_peak1} kDa is the unbound Protein B dimer.")
print(f"Peak 2 at {exp2_peak2} kDa corresponds to a complex of the remaining proteins.")
print(f"Mass Calculation: {A_dimer_mass} (A-dimer) + {C_monomer_mass} (C) + {D_monomer_mass} (D) = {complex_ACD_mass} kDa.")
print("Conclusion 2: Non-phosphorylated Protein A (as a dimer) has a higher affinity for Proteins C and D than the Protein B dimer does.\n")

# 3. Analyze Experiment 3
print("--- Analysis of Experiment 3 (Mixture with Kinase) ---")
exp3_peak1 = 25
exp3_peak3 = 460
complex_BCD_mass = B_dimer_mass + C_monomer_mass + D_monomer_mass
print(f"Peak 1 at {exp3_peak1} kDa is monomeric Protein A. This means phosphorylation breaks the A-dimer and causes it to dissociate from the complex.")
print(f"Peak 3 at {exp3_peak3} kDa corresponds to a new complex.")
print(f"Mass Calculation: {B_dimer_mass} (B-dimer) + {C_monomer_mass} (C) + {D_monomer_mass} (D) = {complex_BCD_mass} kDa.")
print("Conclusion 3: Phosphorylation of Protein A decreases its affinity for C and D, allowing the Protein B dimer to form a complex with them.\n")

# 4. Analyze Experiment 4
print("--- Analysis of Experiment 4 (Dephosphorylation of A) ---")
exp4_peak1 = 50
exp4_peak2 = 460
print(f"Peak 1 at {exp4_peak1} kDa shows that dephosphorylated Protein A re-forms its dimer.")
print(f"Peak 2 at {exp4_peak2} kDa shows the (B-dimer)-C-D complex remains stable.")
print("Conclusion 4: The (B-dimer)-C-D complex is stable once formed.\n")

# 5. Evaluate Answer Choices
print("--- Final Evaluation of Answer Choices ---")
print("A, B, C, D, E, F, H are incorrect based on the conclusions above.")
print("G is the only statement fully supported by the data:")
print("  - 'Protein B never has a higher affinity for protein C or D than nonphosphorylated protein A' (Correct, from Exp 2).")
print("  - 'protein B always exists as a dimer' (Correct, from all experiments).")
print("  - 'in opposition to proteins A, C, and D' (Correct, as A can be a monomer, and C/D are monomers).")

# Final Answer
print("\nFinal Answer is G.")

# Restore original stdout
sys.stdout = original_stdout
# Get the content from the buffer
output_content = output_buffer.getvalue()

# Print the captured output
print(output_content)

print("<<<G>>>")