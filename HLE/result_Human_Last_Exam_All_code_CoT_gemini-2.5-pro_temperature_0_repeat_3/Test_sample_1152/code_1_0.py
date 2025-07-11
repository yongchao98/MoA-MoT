import sys
import io

# Redirect stdout to capture prints
old_stdout = sys.stdout
sys.stdout = captured_output = io.StringIO()

# --- Data Representation ---
# Mortality rates are stored as decimals (e.g., 35% -> 0.35)
control_data = {
    "non_infected_mortality": 0.10
}

exp1_fungus_A = {
    "buck": {"mortality": 0.35},
    "sunflower": {"mortality": 0.10},
    "lavender": {"mortality": 0.35},
    "canola": {"mortality": 0.35},
    "milkweed": {"mortality": 0.35},
    "aster": {"mortality": 0.35},
    "mixed": {"mortality": 0.35},
}

exp2_productivity_A = {
    "buck": {"not_infected": 45, "infected": 10},
    "sunflower": {"not_infected": 30, "infected": 20},
    "lavender": {"not_infected": 30, "infected": 10},
    "canola": {"not_infected": 30, "infected": 8},
    "milkweed": {"not_infected": 30, "infected": 9},
    "aster": {"not_infected": 30, "infected": 11},
    "mixed": {"not_infected": 32, "infected": 12},
}

exp3_fungus_B = {
    "buck": {"mortality": 0.20},
    "sunflower": {"mortality": 0.20},
    "lavender": {"mortality": 0.20},
    "canola": {"mortality": 0.20},
    "milkweed": {"mortality": 0.20},
    "aster": {"mortality": 0.20},
    "mixed": {"mortality": 0.20},
}

exp4_fungus_C = {
    "buck": {"mortality": 0.10},
    "sunflower": {"mortality": 0.10},
    "lavender": {"mortality": 0.10},
    "canola": {"mortality": 0.10},
    "milkweed": {"mortality": 0.10},
    "aster": {"mortality": 0.10},
    "mixed": {"mortality": 0.10},
}

exp5_productivity_C = {
    "buck": {"infected": 60},
    "sunflower": {"infected": 25},
    "lavender": {"infected": 50},
    "canola": {"infected": 50},
    "milkweed": {"infected": 50},
    "aster": {"infected": 50},
    "mixed": {"infected": 52},
}

# --- Step 1 & 2: Analyze Pathogenicity and Commensalism ---
print("--- Analysis of Fungal Types ---")
control_mortality = control_data["non_infected_mortality"]
print(f"The baseline mortality rate for non-infected bees is {control_mortality*100}%.")

# Fungus A
max_mortality_A = max(pollen["mortality"] for pollen in exp1_fungus_A.values())
print(f"Fungus A: Maximum mortality is {max_mortality_A*100}%. Since {max_mortality_A*100}% > {control_mortality*100}%, Fungus A is a pathogen.")

# Fungus B
max_mortality_B = max(pollen["mortality"] for pollen in exp3_fungus_B.values())
print(f"Fungus B: Maximum mortality is {max_mortality_B*100}%. Since {max_mortality_B*100}% > {control_mortality*100}%, Fungus B is a pathogen.")

# Fungus C
max_mortality_C = max(pollen["mortality"] for pollen in exp4_fungus_C.values())
print(f"Fungus C: Maximum mortality is {max_mortality_C*100}%. Since this is equal to the control rate of {control_mortality*100}%, Fungus C does not increase mortality.")

# Check Fungus C productivity
print("Checking productivity for Fungus C infected bees vs non-infected:")
for pollen, data in exp5_productivity_C.items():
    not_infected_eggs = exp2_productivity_A[pollen]["not_infected"]
    infected_eggs = data["infected"]
    if infected_eggs >= not_infected_eggs:
        print(f"- {pollen.capitalize()}: Productivity increased or stayed same ({not_infected_eggs} to {infected_eggs} eggs).")
    else:
        print(f"- {pollen.capitalize()}: Productivity decreased ({not_infected_eggs} to {infected_eggs} eggs).")
print("Since Fungus C does not increase mortality and generally increases productivity, it is a commensal (or mutualist), not a pathogen.")

print("\nCONCLUSION 1: Fungus A and B are pathogens. Fungus C is a commensal.")

# --- Step 3, 4, 5: Evaluate other claims ---
print("\n--- Evaluating Other Claims ---")

# Claim: Buck pollen is always best for productivity
print("\nEvaluating claim: 'Buck pollen always leads to highest productivity'")
infected_A_eggs_buck = exp2_productivity_A["buck"]["infected"]
infected_A_eggs_sunflower = exp2_productivity_A["sunflower"]["infected"]
print(f"When infected with Fungus A, Buck leads to {infected_A_eggs_buck} eggs, while Sunflower leads to {infected_A_eggs_sunflower} eggs.")
print(f"Since {infected_A_eggs_buck} < {infected_A_eggs_sunflower}, the claim that Buck *always* produces more eggs is FALSE.")

# --- Final Conclusion ---
print("\n--- Final Answer Determination ---")
print("The analysis shows that 'Fungus A and B are pathogens. Fungus C is a commensal.' is a correct statement derived from the data.")
print("Other options contain incorrect claims, such as those about pollen productivity being 'always' the best, or incorrect comparisons of deadliness.")
print("Therefore, option I is the most accurate and complete answer.")

final_answer = "I"

# Restore stdout and print the captured output
sys.stdout = old_stdout
print(captured_output.getvalue())
print(f"<<<{final_answer}>>>")