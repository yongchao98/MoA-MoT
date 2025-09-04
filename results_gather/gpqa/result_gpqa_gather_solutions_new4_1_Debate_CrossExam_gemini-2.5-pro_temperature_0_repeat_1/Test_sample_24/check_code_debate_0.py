import sys
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def install_rdkit():
    """Installs the rdkit library if it's not already installed."""
    try:
        import rdkit
    except ImportError:
        print("RDKit library not found. Attempting to install...")
        import subprocess
        subprocess.check_call([sys.executable, "-m", "pip", "install", "rdkit-pypi"])
        print("RDKit installed successfully.")

def get_mol_formula(mol):
    """Calculates the molecular formula from an RDKit molecule object."""
    if mol is None:
        return None
    return rdMolDescriptors.CalcMolFormula(mol)

def check_functional_group(mol, smarts_pattern):
    """Checks if a molecule contains a specific functional group using SMARTS."""
    if mol is None:
        return False
    pattern = Chem.MolFromSmarts(smarts_pattern)
    return mol.HasSubstructMatch(pattern)

def check_answer():
    """
    Checks the correctness of the answer by verifying the chemical logic
    of the transformations using rdkit.
    """
    # Define SMILES strings for all relevant molecules
    # Note: A representative structure is used for the complex diol A,
    # preserving its molecular formula (C12H22O2) and diol nature for logical checks.
    molecules = {
        "product_A": "CC1CCC2(C)CCCC(=O)C2C1",  # 2,8-dimethylspiro[4.5]decan-6-one (C12H20O)
        "product_B": "CC(C)=CC[CH](O)c1ccccc1", # 4-methyl-1-phenylpent-3-en-1-ol (C12H16O)
        "reactant_A_correct": "CC1(O)C(O)C2CCCCC2C1(C)C", # Proxy for 2,7-dimethyloctahydronaphthalene-4a,8a-diol (C12H22O2)
        "reactant_B_correct": "c1ccccc1COOCC=C(C)C", # (((3-methylbut-2-en-1-yl)oxy)methyl)benzene (C12H16O)
        "reactant_A_alt": "CC1CCC2(C)CCCC(O)C2C1", # 2,8-dimethylspiro[4.5]decan-6-ol (C12H22O)
        "reactant_B_alt": "CC(C)=CCC(=O)c1ccccc1", # 4-methyl-1-phenylpent-3-en-1-one (C12H14O)
    }

    # Define the options from the question
    options = {
        "A": {"A": "reactant_A_correct", "B": "reactant_B_correct"},
        "B": {"A": "reactant_A_alt", "B": "reactant_B_alt"}, # Note: This is a mix from the original options for checking purposes
        "C": {"A": "reactant_A_alt", "B": "reactant_B_correct"},
        "D": {"A": "reactant_A_correct", "B": "reactant_B_alt"},
    }

    # The final answer provided by the LLM
    llm_answer = "A"

    # --- Verification Logic ---
    
    # Get the proposed reactants from the LLM's answer
    proposed_option = options.get(llm_answer)
    if not proposed_option:
        return f"Invalid answer option '{llm_answer}'. Please choose from A, B, C, D."

    proposed_reactant_A_key = proposed_option["A"]
    proposed_reactant_B_key = proposed_option["B"]

    # Create RDKit molecule objects
    mols = {key: Chem.MolFromSmiles(smi) for key, smi in molecules.items()}

    # --- Check Reaction A ---
    # Principle: Pinacol rearrangement (Diol -> Ketone + H2O)
    
    # 1. Check the proposed reactant A
    reactant_A_mol = mols[proposed_reactant_A_key]
    product_A_mol = mols["product_A"]
    
    formula_reactant_A = get_mol_formula(reactant_A_mol)
    formula_product_A = get_mol_formula(product_A_mol)
    
    # Expected formula for reactant A is C12H22O2, product is C12H20O
    if formula_reactant_A != "C12H22O2" or formula_product_A != "C12H20O":
        return f"Incorrect molecular formula calculation for Reaction A. Reactant: {formula_reactant_A}, Product: {formula_product_A}"

    # Check if it's a diol (has 2 OH groups)
    if reactant_A_mol.GetSubstructMatches(Chem.MolFromSmarts("[#8H1]")) != 2:
         return f"Reason: Reactant A from option {llm_answer} is not a diol, which is required for a Pinacol rearrangement."

    # 2. Check the alternative reactant A to confirm it's wrong
    alt_reactant_A_mol = mols["reactant_A_alt"]
    if alt_reactant_A_mol.GetSubstructMatches(Chem.MolFromSmarts("[#8H1]")) == 2:
        return "Logic Error: The alternative reactant for A should not be a diol."
    
    # --- Check Reaction B ---
    # Principle: Wittig rearrangement (Ether -> Alcohol, Isomerization)

    # 1. Check the proposed reactant B
    reactant_B_mol = mols[proposed_reactant_B_key]
    product_B_mol = mols["product_B"]

    formula_reactant_B = get_mol_formula(reactant_B_mol)
    formula_product_B = get_mol_formula(product_B_mol)

    # Check for isomerization (same molecular formula)
    if formula_reactant_B != formula_product_B:
        return f"Reason: Reactant B from option {llm_answer} ({formula_reactant_B}) is not an isomer of the product ({formula_product_B}). A Wittig rearrangement is an isomerization."

    # Check functional groups (Ether -> Alcohol)
    if not check_functional_group(reactant_B_mol, "[#6]-[#8]-[#6]"):
        return f"Reason: Reactant B from option {llm_answer} is not an ether, which is required for a Wittig rearrangement."
    if not check_functional_group(product_B_mol, "[#6][#8H1]"):
        return f"Logic Error: Product B should be an alcohol."

    # 2. Check the alternative reactant B to confirm it's a less likely choice
    alt_reactant_B_mol = mols["reactant_B_alt"]
    formula_alt_B = get_mol_formula(alt_reactant_B_mol)
    if formula_alt_B == formula_product_B:
        return "Logic Error: The alternative reactant for B should not be an isomer of the product."
    if not check_functional_group(alt_reactant_B_mol, "[#6][C](=[O])[#6]"):
        return "Logic Error: The alternative reactant for B should be a ketone."

    # If all checks pass, the logic of the answer is sound.
    return "Correct"

# Run the check
try:
    install_rdkit()
    result = check_answer()
    print(result)
except Exception as e:
    print(f"An error occurred during the check: {e}")
