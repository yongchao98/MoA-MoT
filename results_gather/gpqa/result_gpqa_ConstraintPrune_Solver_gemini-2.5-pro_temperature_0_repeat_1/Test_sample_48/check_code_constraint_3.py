# This code block verifies the constraints used in the analysis.

def get_formula(name):
    """A simple dictionary to store molecular formulas for the compounds in the question."""
    formulas = {
        # Reactants
        "1,1-dimethoxyethan-1-amine": "C4H11NO2",
        "but-3-en-2-ol": "C4H8O",
        "(3R,4S)-3,4-dimethylhexa-1,5-diyne": "C8H10",
        "2-((vinyloxy)methyl)but-1-ene": "C7H12O",
        # Products
        "(1Z,2E)-1,2-diethylidenecyclobutane": "C8H12",
        "(3Z,4E)-3,4-diethylidenecyclobut-1-ene": "C8H10",
        "4-methylenehexanal": "C7H12O",
        "4-methylenehexan-1-ol": "C7H14O",
    }
    return formulas.get(name, "Unknown")

def check_functional_group(name):
    """Checks if a name corresponds to a carbonyl (aldehyde/ketone) or alcohol."""
    if name.endswith("al") or name.endswith("one"):
        return "Carbonyl"
    if name.endswith("ol"):
        return "Alcohol"
    return "Other"

print("--- Constraint Check: Reaction B (Isomerization) ---")
reactant_b_formula = get_formula("(3R,4S)-3,4-dimethylhexa-1,5-diyne")
candidate_b1 = "(1Z,2E)-1,2-diethylidenecyclobutane"
candidate_b2 = "(3Z,4E)-3,4-diethylidenecyclobut-1-ene"
print(f"Reactant B is {reactant_b_formula}. Candidate '{candidate_b1}' ({get_formula(candidate_b1)}): Pass/Fail -> {get_formula(candidate_b1) == reactant_b_formula}")
print(f"Reactant B is {reactant_b_formula}. Candidate '{candidate_b2}' ({get_formula(candidate_b2)}): Pass/Fail -> {get_formula(candidate_b2) == reactant_b_formula}")

print("\n--- Constraint Check: Reaction C (Isomerization & Functional Group) ---")
reactant_c_formula = get_formula("2-((vinyloxy)methyl)but-1-ene")
candidate_c1 = "4-methylenehexanal"
candidate_c2 = "4-methylenehexan-1-ol"
# Check candidate 1
formula_pass_c1 = get_formula(candidate_c1) == reactant_c_formula
type_pass_c1 = check_functional_group(candidate_c1) == "Carbonyl"
print(f"Reactant C is {reactant_c_formula}. Candidate '{candidate_c1}' ({get_formula(candidate_c1)}, {check_functional_group(candidate_c1)}): Pass/Fail -> {formula_pass_c1 and type_pass_c1}")
# Check candidate 2
formula_pass_c2 = get_formula(candidate_c2) == reactant_c_formula
type_pass_c2 = check_functional_group(candidate_c2) == "Carbonyl"
print(f"Reactant C is {reactant_c_formula}. Candidate '{candidate_c2}' ({get_formula(candidate_c2)}, {check_functional_group(candidate_c2)}): Pass/Fail -> {formula_pass_c2 and type_pass_c2}")