import sys
import io

# Redirect stdout to capture print statements for final output formatting
old_stdout = sys.stdout
sys.stdout = captured_output = io.StringIO()

# Step 1: Define the data from the problem description.
masses = {
    "GIKS3": 60,
    "CaPK1": 45,
    "CaPK2": 102,
    "CaPK3": 65,
    "CaPK4": 55,
    "CaPK5": 39,
}

# SEC-MALS results (Experiment 1)
sec_mals_results = {
    "GIKS3 alone": [120],
    "GIKS3+CaPK1": [45, 120],
    "GIKS3+CaPK2": [222],
    "GIKS3+CaPK3": [65, 120, 185],
    "GIKS3+CaPK4": [55, 120],
    "GIKS3+CaPK5": [39, 120, 159],
}

# Phosphorylation results (Experiment 2)
phosphorylation_results = {
    "GIKS3-wt+CaPK1": [45],
    "GIKS3-S25A+CaPK1": [45],
    "GIKS3-wt+CaPK2": [60, 102],
    "GIKS3-S25A+CaPK2": [60, 102],
    "GIKS3-wt+CaPK3": [60, 65],
    "GIKS3-S25A+CaPK3": [65],
    "GIKS3-wt+CaPK4": [55, 60],
    "GIKS3-S25A+CaPK4": [55],
    "GIKS3-wt+CaPK5": [],
    "GIS3-S25A+CaPK5": [],
}

# Activity results (Experiment 3)
activity_results = {
    "GIKS3-wt+CaPK1": 0,
    "GIKS3-S25A+CaPK1": 0,
    "GIKS3-wt+CaPK2": 3,
    "GIKS3-S25A+CaPK2": 0,
    "GIKS3-wt+CaPK3": 3,
    "GIKS3-S25A+CaPK3": 0,
    "GIKS3-wt+CaPK4": 3,
    "GIKS3-S25A+CaPK4": 0,
    "GIKS3-wt+CaPK5": 0,
    "GIKS3-S25A+CaPK5": 0,
}

# Dictionary to store derived conclusions
conclusions = {}

print("Step-by-step Analysis of Experimental Data:")
print("=" * 50)

# Step 2: Analysis of Experiment 1: SEC-MALS (Interaction)
print("Analysis of Experiment 1: SEC-MALS (Interaction)")
giks3_mass = masses["GIKS3"]
giks3_dimer_mass = sec_mals_results["GIKS3 alone"][0]
print(f"GIKS3 has a monomer mass of {giks3_mass} kDa. The control peak at {giks3_dimer_mass} kDa indicates GIKS3 forms a dimer.")
print("-" * 20)

for i in range(1, 6):
    kinase = f"CaPK{i}"
    kinase_mass = masses[kinase]
    # The "final equation" for complex formation
    expected_complex_mass = giks3_dimer_mass + kinase_mass
    observed_peaks = sec_mals_results[f"GIKS3+{kinase}"]
    
    interacts = expected_complex_mass in observed_peaks
    
    conclusions[kinase] = {'interacts': interacts}
    
    print(f"Analysis for {kinase}:")
    print(f"  Equation: GIKS3 dimer ({giks3_dimer_mass} kDa) + {kinase} ({kinase_mass} kDa) = Expected Complex ({expected_complex_mass} kDa)")
    print(f"  Observed peaks for GIKS3+{kinase}: {observed_peaks} kDa.")
    if interacts:
        print(f"  Conclusion: {kinase} INTERACTS with GIKS3.")
    else:
        print(f"  Conclusion: {kinase} does NOT form a stable complex with GIKS3 detectable by SEC-MALS.")
    print("-" * 20)

print("\n" + "=" * 50)
# Step 3: Analysis of Experiments 2 & 3: Phosphorylation and Activity
print("Analysis of Experiments 2 & 3: Phosphorylation and Activity")
print("Activation of GIKS3 requires phosphorylation at Serine 25. This is tested by comparing activity with GIKS3-wt vs. GIKS3-S25A.")
print("-" * 20)

for i in range(1, 6):
    kinase = f"CaPK{i}"
    
    # From Experiment 3 (Activity)
    activates_wt = activity_results[f"GIKS3-wt+{kinase}"] > 0
    activates_s25a = activity_results[f"GIKS3-S25A+{kinase}"] > 0
    
    phosphorylates_s25 = activates_wt and not activates_s25a
    conclusions[kinase]['phosphorylates_s25'] = phosphorylates_s25
    
    print(f"Analysis for {kinase}:")
    print(f"  Activity Assay: Rate with GIKS3-wt = {activity_results[f'GIKS3-wt+{kinase}']} mmol/min, Rate with GIKS3-S25A = {activity_results[f'GIKS3-S25A+{kinase}']} mmol/min.")
    if phosphorylates_s25:
        print(f"  Conclusion: {kinase} ACTIVATES GIKS3 by phosphorylating Serine 25.")
    else:
        print(f"  Conclusion: {kinase} does NOT activate GIKS3 via Serine 25.")
    print("-" * 20)

print("\n" + "=" * 50)
# Step 4: Synthesize Conclusions
print("Final Summary of Conclusions:")
for i in range(1, 6):
    kinase = f"CaPK{i}"
    interacts_str = "Yes" if conclusions[kinase]['interacts'] else "No"
    phos_s25_str = "Yes" if conclusions[kinase]['phosphorylates_s25'] else "No"
    print(f"{kinase}: Forms stable complex? {interacts_str}. Activates via S25? {phos_s25_str}.")

print("=" * 50)
# Step 5: Evaluate Answer Choices
print("\nEvaluating Answer Choices:")

# A. CaPK2 and CaPK3 can phosphorylate GIKS3 on serine 25. CaPK4 does not interact with GIKS3 and CaPK1 does not interact with GIKS3.
a_part1 = conclusions['CaPK2']['phosphorylates_s25'] and conclusions['CaPK3']['phosphorylates_s25']
a_part2 = not conclusions['CaPK4']['interacts']
a_part3 = not conclusions['CaPK1']['interacts']
is_a_correct = a_part1 and a_part2 and a_part3
print(f"Choice A: 'CaPK2 and CaPK3 can phosphorylate GIKS3 on S25' (True). 'CaPK4 does not interact' (True). 'CaPK1 does not interact' (True). -> Overall: {is_a_correct}")

# B. Only CaPK3 and CaPK4 can activate GIKS3. The complex between CaPK4 and GIKS3 was not detected in the SEC-MALS experiment.
b_part1 = (not conclusions['CaPK2']['phosphorylates_s25']) and conclusions['CaPK3']['phosphorylates_s25'] and conclusions['CaPK4']['phosphorylates_s25']
is_b_correct = b_part1
print(f"Choice B: 'Only CaPK3 and CaPK4 can activate' (False, because CaPK2 also activates). -> Overall: {is_b_correct}")

# I. Only CaPK2 and CaPK3 can phosphorylate GIKS3 on serine 25.
i_part1 = conclusions['CaPK2']['phosphorylates_s25'] and conclusions['CaPK3']['phosphorylates_s25'] and (not conclusions['CaPK4']['phosphorylates_s25'])
is_i_correct = i_part1
print(f"Choice I: 'Only CaPK2 and CaPK3 can phosphorylate' (False, because CaPK4 also does). -> Overall: {is_i_correct}")

print("\n" + "=" * 50)
print("Final Decision:")
print("Choice A is composed of statements that are all factually correct based on the provided data.")
print("Other choices contain inaccuracies (e.g., using 'Only' when multiple kinases are involved) or are less complete summaries.")
print("The statement highlights key findings: successful phosphorylation by CaPK2 and CaPK3, and the lack of stable interaction for CaPK1 and CaPK4 (a particularly interesting finding for the active CaPK4).")

final_answer = "A"
print(f"\nThe most accurate and complete conclusion is provided in choice A.")
print(f"<<<{final_answer}>>>")

# Restore stdout and print the captured output
sys.stdout = old_stdout
print(captured_output.getvalue())