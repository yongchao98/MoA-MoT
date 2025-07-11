# Define protein masses based on sequence (theoretical) and experiment 1 (observed oligomeric state)
theo_mass = {'A': 25, 'B': 150, 'C': 60, 'D': 100}
exp1_mass = {'A': 50, 'B': 300, 'C': 60, 'D': 100}
kinase_mass = 40

print("--- Step 1: Analyze Oligomeric States from Experiment 1 ---")
# Calculate oligomeric states
state_A = exp1_mass['A'] / theo_mass['A']
state_B = exp1_mass['B'] / theo_mass['B']
state_C = exp1_mass['C'] / theo_mass['C']
state_D = exp1_mass['D'] / theo_mass['D']
print(f"Protein A (native): {exp1_mass['A']} kDa / {theo_mass['A']} kDa = {int(state_A)}. It is a homodimer.")
print(f"Protein B (native): {exp1_mass['B']} kDa / {theo_mass['B']} kDa = {int(state_B)}. It is a homodimer.")
print(f"Protein C (native): {exp1_mass['C']} kDa / {theo_mass['C']} kDa = {int(state_C)}. It is a monomer.")
print(f"Protein D (native): {exp1_mass['D']} kDa / {theo_mass['D']} kDa = {int(state_D)}. It is a monomer.\n")

print("--- Step 2: Analyze Complex Formation in Experiment 2 ---")
# Experiment 2 results: one peak at 300 kDa, one at 210 kDa
exp2_peak1 = 300
exp2_peak2 = 210
# Identify the components of the 210 kDa peak
complex_ACD_mass = exp1_mass['A'] + exp1_mass['C'] + exp1_mass['D']
print(f"Observed peaks are {exp2_peak1} kDa and {exp2_peak2} kDa.")
print(f"The {exp2_peak1} kDa peak corresponds to unbound Protein B dimer.")
print(f"Let's calculate the mass of a complex of the remaining proteins (A, C, D):")
print(f"Protein A (dimer) + Protein C + Protein D = {exp1_mass['A']} kDa + {exp1_mass['C']} kDa + {exp1_mass['D']} kDa = {complex_ACD_mass} kDa.")
print("This matches the second peak. Conclusion: Non-phosphorylated Protein A outcompetes Protein B for binding to C and D.\n")

print("--- Step 3: Analyze Effect of Phosphorylation in Experiment 3 ---")
# Experiment 3 results: peaks at 25, 40, and 460 kDa
exp3_peak1 = 25
exp3_peak2 = 40
exp3_peak3 = 460
# Identify the 460 kDa complex
complex_BCD_mass = exp1_mass['B'] + exp1_mass['C'] + exp1_mass['D']
print(f"Observed peaks are {exp3_peak1} kDa, {exp3_peak2} kDa (Kinase), and {exp3_peak3} kDa.")
print(f"The {exp3_peak1} kDa peak is monomeric Protein A, indicating phosphorylation caused its dimer to dissociate.")
print("Let's calculate the mass of a complex of the remaining proteins (B, C, D):")
print(f"Protein B (dimer) + Protein C + Protein D = {exp1_mass['B']} kDa + {exp1_mass['C']} kDa + {exp1_mass['D']} kDa = {complex_BCD_mass} kDa.")
print("This matches the third peak. Conclusion: Phosphorylation of A removes it from competition, allowing B, C, and D to form a complex.\n")

print("--- Step 4: Final Conclusion ---")
print("Based on the experiments:")
print("1. (From Exp 2) Non-phosphorylated Protein A has a higher affinity for Proteins C and D than Protein B does.")
print("2. (From Exp 1, 2, 3, 4) Protein B is consistently a homodimer (300 kDa).")
print("3. (From Exp 1, 3) Protein A can be a monomer or a dimer, while C and D are monomers. So, Protein B is the only one that is *always* a dimer.")
print("These points directly support answer choice G.\n")
print("Evaluating Answer G: 'Protein B never has a higher affinity for protein C or D than nonphosphorylated protein A, protein B always exists as a dimer in opposition to proteins A, C, and D.'")
print("This statement is fully consistent with our analysis.")

<<<G>>>