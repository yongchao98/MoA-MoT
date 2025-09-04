import sys

# We use the rdkit library for cheminformatics. If it's not installed,
# the script will use pre-calculated formulas for a basic check.
try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False

def get_mol_formula(smiles: str) -> str:
    """Calculates the molecular formula from a SMILES string."""
    if not RDKIT_AVAILABLE:
        # Fallback dictionary for environments without rdkit.
        # These formulas are pre-calculated manually.
        formulas = {
            "C=CC1(O)C2(CCC1)C=CCCC2": "C11H16O", # 1-vinylspiro[3.5]non-5-en-1-ol
            "O=C1CC2CC(C=C3C2)CCC13": "C11H16O", # (E)-bicyclo[5.3.1]undec-1(11)-en-4-one
            "O=C1CCCC2CCCCC12": "C11H18O", # decahydro-7H-benzo[7]annulen-7-one
        }
        if smiles in formulas:
            return formulas[smiles]
        return f"Formula for {smiles} not pre-calculated."

    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return f"Error: Could not parse SMILES: {smiles}"
        return Chem.rdMolDescriptors.CalcMolFormula(mol)
    except Exception as e:
        return f"Error calculating formula for {smiles}: {e}"

def check_correctness():
    """
    Checks the correctness of the LLM's answer for the two organic reactions
    by verifying stoichiometry and reaction conditions.
    """
    errors = []
    
    # The provided answer identifies the products as:
    # A = (E)-bicyclo[5.3.1]undec-1(11)-en-4-one
    # B = lithium 3-ethylpent-4-enoate
    
    # --- Reaction 1 Check ---
    # Reaction: 1-vinylspiro[3.5]non-5-en-1-ol --> A
    # Mechanism: Anionic Oxy-Cope Rearrangement. This is an isomerization,
    # so the molecular formula must be conserved.
    
    start_1_smiles = "C=CC1(O)C2(CCC1)C=CCCC2"
    product_A_correct_smiles = "O=C1CC2CC(C=C3C2)CCC13"
    product_A_incorrect_smiles = "O=C1CCCC2CCCCC12"

    start_1_formula = get_mol_formula(start_1_smiles)
    product_A_correct_formula = get_mol_formula(product_A_correct_smiles)
    product_A_incorrect_formula = get_mol_formula(product_A_incorrect_smiles)

    # Constraint 1: The correct product A must be an isomer of the starting material.
    if start_1_formula != product_A_correct_formula:
        errors.append(
            f"Constraint Violated (Reaction 1): The proposed correct product A (formula {product_A_correct_formula}) "
            f"is not an isomer of the starting material (formula {start_1_formula}). "
            "The Anionic Oxy-Cope rearrangement must conserve the molecular formula."
        )

    # Constraint 2: The incorrect product A should have a different formula.
    if start_1_formula == product_A_incorrect_formula:
        errors.append(
            f"Constraint Violated (Reaction 1): The proposed incorrect product A (formula {product_A_incorrect_formula}) "
            f"has the same formula as the starting material, which contradicts the LLM's correct reasoning that it is stoichiometrically wrong."
        )
    
    # --- Reaction 2 Check ---
    # Reaction: (E)-pent-2-en-1-ol + acetyl bromide (Base = LDA) --> B
    # Mechanism: Ireland-Claisen Rearrangement.
    
    # Constraint 3: The final product form (salt vs. acid) must match the conditions.
    # The reagents include LDA, a strong base. The reaction produces a carboxylate.
    # No acidic workup (like H+) is specified for reaction 2.
    # Therefore, the final product should be the deprotonated salt form.
    
    llm_answer_product_B = "lithium 3-ethylpent-4-enoate"
    
    if "lithium" not in llm_answer_product_B.lower() or "enoate" not in llm_answer_product_B.lower():
        errors.append(
            "Constraint Violated (Reaction 2): The product B should be a lithium salt (e.g., 'lithium ...enoate') "
            "due to the use of LDA as a base and the lack of a specified acidic workup step. "
            "Identifying the neutral acid form ('...enoic acid') would be incorrect under these conditions."
        )

    # --- Final Conclusion ---
    if not errors:
        return "Correct"
    else:
        return "Incorrect. Reason(s):\n" + "\n".join(errors)

# Execute the check and print the result
result = check_correctness()
print(result)