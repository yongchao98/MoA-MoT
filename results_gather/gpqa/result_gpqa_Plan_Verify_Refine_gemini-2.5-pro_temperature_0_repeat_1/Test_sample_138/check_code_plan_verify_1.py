from rdkit import Chem
from rdkit.Chem import AllChem

def get_canonical_smiles(name_or_smiles):
    """Converts a chemical name or SMILES string to a canonical SMILES string."""
    try:
        # A simple dictionary to map names to SMILES strings
        name_to_smiles = {
            # Reactants
            "4-isopropylcyclohexan-1-one": "CC(C)C1CCC(=O)CC1",
            "5-methylhexane-2,3-diol": "CCC(C)C(O)C(O)C",
            "4-isopropyl-2-methoxycyclohexan-1-ol": "COC1C(C(C)C)CCC(O)C1",
            "5-methylhexan-2-one": "CC(C)CCC(=O)C",
            # Products
            "4-isopropylcyclohexane-1,2-dione": "CC(C)C1CCC(=O)C(=O)C1",
            "5-methylhexane-2,3-dione": "CC(C)CC(=O)C(=O)C"
        }
        smiles = name_to_smiles.get(name_or_smiles, name_or_smiles)
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        return Chem.MolToSmiles(mol, canonical=True)
    except:
        return None

def check_reaction(reactant_name, product_name):
    """
    Checks if a reactant can form a product via alpha-oxidation of a ketone.
    Returns (True, "") if correct, or (False, "reason") if incorrect.
    """
    reactant_smiles = get_canonical_smiles(reactant_name)
    product_smiles = get_canonical_smiles(product_name)
    
    reactant_mol = Chem.MolFromSmiles(reactant_smiles)

    # Constraint 1: The reactant must be a ketone.
    ketone_pattern = Chem.MolFromSmarts('[#6][C](=O)[#6]')
    if not reactant_mol.HasSubstructMatch(ketone_pattern):
        return False, f"Reactant '{reactant_name}' is not a ketone, which is required for this reaction."

    # Constraint 2: The ketone must have an alpha-methylene that can be oxidized.
    # We will generate all possible alpha-oxidation products and see if the target product is among them.
    possible_products = set()
    # Define the reaction: [CH2:1][C:2](=O)>>[C:1](=O)[C:2](=O)
    rxn = AllChem.ReactionFromSmarts('[CH2:1][C:2](=O)>>[C:1](=O)[C:2](=O)')
    
    product_sets = rxn.RunReactants((reactant_mol,))
    
    if not product_sets:
        return False, f"Reactant '{reactant_name}' does not have an alpha-methylene group to oxidize."

    for products in product_sets:
        for p in products:
            try:
                Chem.SanitizeMol(p)
                possible_products.add(Chem.MolToSmiles(p, canonical=True))
            except:
                continue # Ignore invalid chemical structures

    if product_smiles in possible_products:
        return True, ""
    else:
        return False, f"Oxidation of '{reactant_name}' does not yield '{product_name}'. Possible products: {list(possible_products)}"

def check_correctness():
    """
    This function checks the correctness of the LLM's answer.
    """
    llm_answer = "C"
    
    # Define the question's reactions and options
    target_product_A = "4-isopropylcyclohexane-1,2-dione"
    target_product_B = "5-methylhexane-2,3-dione"
    
    options = {
        "A": ("4-isopropylcyclohexan-1-one", "5-methylhexane-2,3-diol"),
        "B": ("4-isopropyl-2-methoxycyclohexan-1-ol", "5-methylhexane-2,3-diol"),
        "C": ("4-isopropylcyclohexan-1-one", "5-methylhexan-2-one"),
        "D": ("4-isopropyl-2-methoxycyclohexan-1-ol", "5-methylhexan-2-one")
    }

    if llm_answer not in options:
        return f"Invalid answer choice '{llm_answer}'. Must be one of {list(options.keys())}."

    # Get the reactants for the chosen answer
    selected_reactant_A, selected_reactant_B = options[llm_answer]

    # Check the reaction for compound A
    is_A_correct, reason_A = check_reaction(selected_reactant_A, target_product_A)
    if not is_A_correct:
        return f"Incorrect. The proposed reactant A ('{selected_reactant_A}') is wrong. Reason: {reason_A}"

    # Check the reaction for compound B
    is_B_correct, reason_B = check_reaction(selected_reactant_B, target_product_B)
    if not is_B_correct:
        return f"Incorrect. The proposed reactant B ('{selected_reactant_B}') is wrong. Reason: {reason_B}"

    # If both checks pass, the answer is correct.
    return "Correct"

# Run the check
result = check_correctness()
print(result)