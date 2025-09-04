from rdkit import Chem
from rdkit.Chem import Descriptors

def check_correctness():
    """
    Checks the correctness of the selected answer for the two Michael addition reactions.
    It verifies mass conservation, regioselectivity, and chemical plausibility of the product structures.
    """
    try:
        # The final answer provided by the LLM analysis.
        final_answer = "A"

        # --- Data Representation ---
        # SMILES strings for reactants
        reactant_A1_smiles = "COC(=O)C1C(=O)CCCC1"  # methyl 2-oxocyclohexane-1-carboxylate
        reactant_A2_smiles = "CC1=CC(=C(C=C1)S(=O)C=C)C"  # 2,4-dimethyl-1-(vinylsulfinyl)benzene

        reactant_B1_smiles = "CCOC(=O)C(CC)CC"  # ethyl 2-ethylbutanoate
        reactant_B2_smiles = "COC(=O)C(=C1CCCC1)C2=CC=CC=C2"  # methyl 2-cyclopentylidene-2-phenylacetate

        # Proposed products from the selected answer (Option A)
        product_A_name = "methyl 1-(2-((2,4-dimethylphenyl)sulfinyl)ethyl)-2-oxocyclohexane-1-carboxylate"
        product_A_smiles = "COC(=O)C1(C(=O)CCCC1)CCS(=O)C2=CC(=C(C=C2)C)C"

        product_B_name = "ethyl 2-ethyl-2-(1-(2-methoxy-2-oxo-1-phenylethyl)cyclopentyl)butanoate"
        product_B_smiles = "CCOC(=O)C(CC)(CC)C1(CCCC1)C(C(=O)OC)C2=CC=CC=C2"

        # --- Verification Checks ---

        # === Check Reaction A ===

        # 1. Mass Conservation for Reaction A
        mol_react_A1 = Chem.MolFromSmiles(reactant_A1_smiles)
        mol_react_A2 = Chem.MolFromSmiles(reactant_A2_smiles)
        mol_prod_A = Chem.MolFromSmiles(product_A_smiles)

        formula_react_A = Descriptors.rdMolDescriptors.CalcMolFormula(Chem.MolFromSmiles(f"{reactant_A1_smiles}.{reactant_A2_smiles}"))
        formula_prod_A = Descriptors.rdMolDescriptors.CalcMolFormula(mol_prod_A)

        if formula_react_A != formula_prod_A:
            return (f"Incorrect: Mass is not conserved in Reaction A. "
                    f"Sum of reactants' formula ({formula_react_A}) does not match product's formula ({formula_prod_A}).")

        # 2. Regioselectivity for Reaction A
        # The name must indicate substitution at C1, not C3.
        if not ("1-(" in product_A_name and "2-oxocyclohexane-1-carboxylate" in product_A_name):
            return (f"Incorrect: Product A's name '{product_A_name}' indicates incorrect regioselectivity. "
                    f"The Michael addition must occur at the most acidic C1 position.")

        # === Check Reaction B ===

        # 3. Mass Conservation for Reaction B
        mol_react_B1 = Chem.MolFromSmiles(reactant_B1_smiles)
        mol_react_B2 = Chem.MolFromSmiles(reactant_B2_smiles)
        mol_prod_B = Chem.MolFromSmiles(product_B_smiles)

        formula_react_B = Descriptors.rdMolDescriptors.CalcMolFormula(Chem.MolFromSmiles(f"{reactant_B1_smiles}.{reactant_B2_smiles}"))
        formula_prod_B = Descriptors.rdMolDescriptors.CalcMolFormula(mol_prod_B)

        if formula_react_B != formula_prod_B:
            return (f"Incorrect: Mass is not conserved in Reaction B. "
                    f"Sum of reactants' formula ({formula_react_B}) does not match product's formula ({formula_prod_B}).")

        # 4. Chemical Plausibility for Reaction B
        # The product of a Michael addition is not a succinate (1,4-dicarbonyl).
        if "succinate" in product_B_name.lower():
            return (f"Incorrect: Product B's name '{product_B_name}' identifies it as a succinate derivative. "
                    f"A Michael addition produces a 1,5-dicarbonyl relationship, not a 1,4-dicarbonyl (succinate). "
                    f"This structure is chemically incorrect for the given reaction.")

        # If all checks pass, the answer is correct.
        return "Correct"

    except Exception as e:
        return f"An error occurred during the check: {e}"

# Execute the check and print the result
result = check_correctness()
print(result)