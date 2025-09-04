import re

def get_chemical_info(name):
    """
    A simplified database to store molecular formula and functional group
    for the chemicals mentioned in the question.
    """
    db = {
        # Products
        "2,8-dimethylspiro[4.5]decan-6-one": {"formula": "C12H20O1", "group": "ketone"},
        "4-methyl-1-phenylpent-3-en-1-ol": {"formula": "C12H16O1", "group": "alcohol"},
        # Potential Reactants
        "2,8-dimethylspiro[4.5]decan-6-ol": {"formula": "C12H22O1", "group": "alcohol"},
        "2,7-dimethyloctahydronaphthalene-4a,8a-diol": {"formula": "C12H22O2", "group": "diol"},
        "4-methyl-1-phenylpent-3-en-1-one": {"formula": "C12H14O1", "group": "ketone"},
        "(((3-methylbut-2-en-1-yl)oxy)methyl)benzene": {"formula": "C12H16O1", "group": "ether"},
    }
    return db.get(name)

def parse_formula(formula_str):
    """Parses a chemical formula string into a dictionary of atom counts."""
    atoms = {}
    # Find all atom-count pairs, e.g., C12, H22, O2
    for match in re.finditer(r'([A-Z][a-z]?)(\d*)', formula_str):
        atom, count = match.groups()
        atoms[atom] = int(count) if count else 1
    return atoms

def check_answer():
    """
    Checks the correctness of the LLM's answer by verifying the chemical logic
    of the proposed reactions.
    """
    llm_answer = "D"
    
    options = {
        "A": {"A": "2,8-dimethylspiro[4.5]decan-6-ol", "B": "4-methyl-1-phenylpent-3-en-1-one"},
        "B": {"A": "2,7-dimethyloctahydronaphthalene-4a,8a-diol", "B": "4-methyl-1-phenylpent-3-en-1-one"},
        "C": {"A": "2,8-dimethylspiro[4.5]decan-6-ol", "B": "(((3-methylbut-2-en-1-yl)oxy)methyl)benzene"},
        "D": {"A": "2,7-dimethyloctahydronaphthalene-4a,8a-diol", "B": "(((3-methylbut-2-en-1-yl)oxy)methyl)benzene"}
    }

    # --- Define Products ---
    product_A_info = get_chemical_info("2,8-dimethylspiro[4.5]decan-6-one")
    product_B_info = get_chemical_info("4-methyl-1-phenylpent-3-en-1-ol")

    # --- Get the proposed reactants from the LLM's answer ---
    chosen_reactants = options.get(llm_answer)
    if not chosen_reactants:
        return f"The answer '{llm_answer}' is not a valid option."

    reactant_A_name = chosen_reactants["A"]
    reactant_B_name = chosen_reactants["B"]
    
    reactant_A_info = get_chemical_info(reactant_A_name)
    reactant_B_info = get_chemical_info(reactant_B_name)

    # --- Verification for Reaction 1 (Pinacol Rearrangement) ---
    # Constraint 1: Reactant must be a diol.
    if reactant_A_info["group"] != "diol":
        return (f"Incorrect. The answer proposes reactant A is '{reactant_A_name}', which is an {reactant_A_info['group']}. "
                f"A Pinacol rearrangement requires a diol as the starting material.")

    # Constraint 2: Stoichiometry must correspond to loss of H2O.
    reactant_atoms = parse_formula(reactant_A_info["formula"])
    product_atoms = parse_formula(product_A_info["formula"])
    
    if not (reactant_atoms.get('C') == product_atoms.get('C') and
            reactant_atoms.get('H') == product_atoms.get('H', 0) + 2 and
            reactant_atoms.get('O') == product_atoms.get('O', 0) + 1):
        return (f"Incorrect. The transformation from reactant A ('{reactant_A_name}', {reactant_A_info['formula']}) "
                f"to the product ({product_A_info['formula']}) does not match the stoichiometry of a Pinacol rearrangement, "
                f"which involves the loss of one water molecule (H2O).")

    # --- Verification for Reaction 2 ([1,2]-Wittig Rearrangement) ---
    # Constraint 1: Reactant must be an ether.
    if reactant_B_info["group"] != "ether":
        # Check the alternative: reduction of a ketone
        if reactant_B_info["group"] == "ketone":
             return (f"Incorrect. The answer proposes reactant B is '{reactant_B_name}', a ketone. "
                     f"While reduction is possible, the primary reaction of a ketone with BuLi (a strong nucleophile) is addition, not reduction. "
                     f"A Wittig rearrangement of an ether is a much more plausible explanation for the specific product formed.")
        return (f"Incorrect. The answer proposes reactant B is '{reactant_B_name}', which is a {reactant_B_info['group']}. "
                f"A [1,2]-Wittig rearrangement requires an ether as the starting material.")

    # Constraint 2: Stoichiometry must be an isomerization (same molecular formula).
    if reactant_B_info["formula"] != product_B_info["formula"]:
        return (f"Incorrect. A Wittig rearrangement is an isomerization, meaning the reactant and product should have the same molecular formula. "
                f"However, the proposed reactant B ('{reactant_B_name}', {reactant_B_info['formula']}) does not match the product's formula ({product_B_info['formula']}).")

    return "Correct"

# Execute the check and print the result
result = check_answer()
print(result)