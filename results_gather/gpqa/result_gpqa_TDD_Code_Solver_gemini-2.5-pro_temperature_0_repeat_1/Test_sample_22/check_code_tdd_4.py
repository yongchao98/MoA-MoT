import sys
from collections import defaultdict

# It's good practice to try importing the necessary library and provide a helpful message if it's not installed.
try:
    from rdkit import Chem
    from rdkit.Chem import rdMolDescriptors
except ImportError:
    print("RDKit library not found.")
    print("Please install it using: pip install rdkit-pypi")
    sys.exit(1)

def get_molecular_formula(smiles: str) -> str:
    """Calculates the molecular formula for a given SMILES string."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return f"Error: Invalid SMILES string '{smiles}'"
    return rdMolDescriptors.CalcMolFormula(mol)

def check_answer():
    """
    Checks the correctness of the given answer by verifying molecular formulas.
    The reaction is an intramolecular cyclization, which is an isomerization.
    Therefore, the products must have the same molecular formula as the reactant.
    """
    # --- Define Reactant and Products using SMILES strings ---

    # Reactant: ((2,2-dimethylbut-3-en-1-yl)oxy)benzene
    reactant_smiles = "c1ccccc1OCC(C)(C)C=C"

    # Options with corresponding product SMILES strings
    options = {
        "A": [
            # (4-bromo-2,2-dimethylbutoxy)benzene - Anti-Markovnikov HBr addition
            "c1ccccc1OCC(C)(C)C(Br)C",
            # (3-bromo-2,2-dimethylbutoxy)benzene - Markovnikov HBr addition
            "c1ccccc1OCC(C)(C)CH(Br)C"
        ],
        "B": [
            # 3,3,4-trimethylchromane - Product of direct cyclization
            "CC1C(C)(C)Cc2c(O1)cccc2",
            # 3-isopropyl-3-methyl-2,3-dihydrobenzofuran - Product of rearrangement and cyclization
            "CC(C)C1(C)Cc2c(O1)cccc2"
        ],
        "C": [
            # (4-bromo-2,2-dimethylbutoxy)benzene
            "c1ccccc1OCC(C)(C)C(Br)C",
            # ((2,3-dimethylbut-2-en-1-yl)oxy)benzene - Isomer of reactant
            "c1ccccc1OCC(C)=C(C)C"
        ],
        "D": [
            # 2-(2,2-dimethylbutyl)phenol - Product of rearrangement (e.g., Claisen) and reduction
            "CCCC(C)(C)c1ccccc1O",
            # 4-(2,2-dimethylbutyl)phenol
            "CCCC(C)(C)c1ccc(O)cc1"
        ]
    }

    # The provided answer from the LLM
    llm_answer = "B"

    # --- Verification Logic ---

    # 1. Calculate the molecular formula of the reactant
    reactant_formula = get_molecular_formula(reactant_smiles)
    
    # 2. Get the products from the chosen answer
    chosen_products_smiles = options.get(llm_answer)
    if not chosen_products_smiles:
        return f"Invalid answer key '{llm_answer}' provided."

    # 3. Check if the products in the chosen answer are isomers of the reactant
    for i, smiles in enumerate(chosen_products_smiles):
        product_formula = get_molecular_formula(smiles)
        if product_formula != reactant_formula:
            return (f"Incorrect. The answer is {llm_answer}, but its proposed products are not all isomers of the reactant.\n"
                    f"Reactant formula: {reactant_formula}\n"
                    f"Product {i+1} ('{smiles}') formula: {product_formula}\n"
                    f"Reason: An intramolecular cyclization is an isomerization reaction, so the products must have the same molecular formula as the reactant.")

    # 4. Verify that other options are incorrect to be thorough
    for option_key, products_smiles in options.items():
        if option_key == llm_answer:
            continue
        
        # Check if any product in the other options incorrectly matches the reactant formula
        # or, more importantly, if they have formulas inconsistent with isomerization.
        is_isomer_option = all(get_molecular_formula(s) == reactant_formula for s in products_smiles)
        if is_isomer_option:
            return (f"Check inconclusive. Option {option_key} also contains products that are isomers of the reactant. "
                    f"A deeper chemical analysis is required beyond atom counting.")

    # If all checks pass, the answer is correct based on this fundamental principle.
    return "Correct"

# Run the check and print the result
result = check_answer()
print(result)