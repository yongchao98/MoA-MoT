import sys

try:
    from rdkit import Chem
    from rdkit.Chem.rdMolDescriptors import CalcMolFormula
except ImportError:
    print("RDKit is not installed. Please install it using 'pip install rdkit-pypi'")
    sys.exit(1)

def check_answer():
    """
    Checks the correctness of the given answer by verifying the chemical transformations.
    """
    # The answer provided by the other LLM
    llm_answer = 'A'

    # --- Data Representation ---
    # We represent the reactants and potential products using SMILES strings.
    # This allows for programmatic analysis of their properties.
    reactants = {
        'B': 'C#CC(C)C(C)C#C',  # (3R,4S)-3,4-dimethylhexa-1,5-diyne
        'C': 'CCC(C=C)COC=C'    # 2-((vinyloxy)methyl)but-1-ene
    }

    # SMILES for the products listed in the options
    product_smiles = {
        'B_option_AC': 'CC=C1C(=CC1)=CC', # (3Z,4E)-3,4-diethylidenecyclobut-1-ene
        'B_option_BD': 'CC=C1CCC(C1)=CC', # (1Z,2E)-1,2-diethylidenecyclobutane
        'C_option_AD': 'CCC(=C)CCC=O',    # 4-methylenehexanal
        'C_option_BC': 'CCC(=C)CCCCO'     # 4-methylenehexan-1-ol
    }

    # Mapping options to their respective product SMILES
    options_map = {
        'A': {'B': product_smiles['B_option_AC'], 'C': product_smiles['C_option_AD']},
        'B': {'B': product_smiles['B_option_BD'], 'C': product_smiles['C_option_BC']},
        'C': {'B': product_smiles['B_option_AC'], 'C': product_smiles['C_option_BC']},
        'D': {'B': product_smiles['B_option_BD'], 'C': product_smiles['C_option_AD']},
    }

    # Get the specific products for the given answer
    answer_products = options_map.get(llm_answer)
    if not answer_products:
        return f"Invalid answer key '{llm_answer}' provided."

    # --- Constraint 1: Check Reaction B (Molecular Formula Conservation) ---
    # Sigmatropic rearrangements must conserve the molecular formula.
    reactant_b_mol = Chem.MolFromSmiles(reactants['B'])
    reactant_b_formula = CalcMolFormula(reactant_b_mol)

    product_b_mol = Chem.MolFromSmiles(answer_products['B'])
    product_b_formula = CalcMolFormula(product_b_mol)

    if reactant_b_formula != product_b_formula:
        return (f"Incorrect. Constraint (Reaction B): Molecular formula is not conserved. "
                f"Reactant formula is {reactant_b_formula}, but product B's formula is {product_b_formula}.")

    # --- Constraint 2: Check Reaction C (Claisen Rearrangement Product) ---
    # A Claisen rearrangement of an allyl vinyl ether must yield a gamma,delta-unsaturated aldehyde/ketone.
    
    # First, check formula conservation as a sanity check for the rearrangement
    reactant_c_mol = Chem.MolFromSmiles(reactants['C'])
    reactant_c_formula = CalcMolFormula(reactant_c_mol)
    
    product_c_mol = Chem.MolFromSmiles(answer_products['C'])
    product_c_formula = CalcMolFormula(product_c_mol)

    if reactant_c_formula != product_c_formula:
        return (f"Incorrect. Constraint (Reaction C): Molecular formula is not conserved. "
                f"Reactant formula is {reactant_c_formula}, but product C's formula is {product_c_formula}.")

    # Second, check for the correct functional group (aldehyde).
    # SMARTS pattern for an aldehyde group: [CX3H1](=O)[#6]
    aldehyde_pattern = Chem.MolFromSmarts('[CX3H1](=O)[#6]')
    if not product_c_mol.HasSubstructMatch(aldehyde_pattern):
        return (f"Incorrect. Constraint (Reaction C): The product of a Claisen rearrangement of an allyl vinyl ether should be an aldehyde. "
                f"The proposed product C is not an aldehyde.")

    # If all checks pass, the answer is correct.
    return "Correct"

# Run the check and print the result
result = check_answer()
print(result)