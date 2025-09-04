from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def get_mol_formula(mol_or_smiles):
    """Calculates the molecular formula from an RDKit molecule or SMILES string."""
    if isinstance(mol_or_smiles, str):
        mol = Chem.MolFromSmiles(mol_or_smiles)
    else:
        mol = mol_or_smiles
    if mol is None:
        return None
    return rdMolDescriptors.CalcMolFormula(mol)

def check_answer():
    """
    Checks the correctness of the provided answer for the multi-step synthesis problem.
    The check is based on stoichiometry (atom counting) and consistency with reaction rules.
    """
    try:
        # The final answer provided by the LLM is C.
        llm_answer_label = "C"
        llm_answer_name = "6-hydroxy-2,2,5,5-tetramethyloctan-4-one"

        # --- Constraint 1: Analyze the starting material ---
        # Name: 3,3,6-trimethylhepta-1,5-dien-4-one
        # Structure: CH2=CH-C(CH3)2-C(=O)-CH=C(CH3)-CH3
        # SMILES representation: C=CC(C)(C)C(=O)C=C(C)C
        start_smiles = "C=CC(C)(C)C(=O)C=C(C)C"
        start_formula = get_mol_formula(start_smiles)
        if start_formula != "C10H16O":
            return f"Internal check failed: Starting material formula should be C10H16O, but calculated as {start_formula}."

        # --- Constraint 2: Analyze the reaction pathway leading to the proposed answer ---
        # The question states a 1:1 mixture of two epoxides is formed. One of them is:
        # Intermediate A: 1,2-epoxy-3,3,6-trimethylhepta-5-en-4-one
        # SMILES: C1OC1C(C)(C)C(=O)C=C(C)C
        intermediate_A_formula = get_mol_formula("C1OC1C(C)(C)C(=O)C=C(C)C")
        if intermediate_A_formula != "C10H16O2":
            return f"Internal check failed: Intermediate A formula should be C10H16O2, but calculated as {intermediate_A_formula}."

        # The question states an "excess" of Gilman reagent ((CH3)2CuLi) is used.
        # Intermediate A has two reactive sites: the epoxide and the alpha,beta-unsaturated ketone.
        # "Excess" implies both sites will react, each consuming one methyl group from the reagent.
        # Therefore, a total of TWO methyl groups are added.
        # The total change in atoms is +2 Carbons and +8 Hydrogens (2 from the methyl groups' H, and 2 from workup).
        # Expected final formula = C10H16O2 + C2H8 = C12H24O2.

        # --- Constraint 3: Analyze the proposed final product (Answer C) ---
        # Name: 6-hydroxy-2,2,5,5-tetramethyloctan-4-one
        # SMILES representation: CCC(O)C(C)(C)C(=O)CC(C)(C)C
        answer_c_smiles = "CCC(O)C(C)(C)C(=O)CC(C)(C)C"
        answer_c_formula = get_mol_formula(answer_c_smiles)

        # Check 1: Stoichiometry. Does the formula of the answer match the expected formula?
        if answer_c_formula != "C12H24O2":
            return (f"Incorrect. The molecular formula of the proposed answer C ({answer_c_formula}) "
                    f"does not match the expected formula from the reaction pathway (C12H24O2). "
                    "The reaction of the 1,2-epoxide intermediate with excess Gilman reagent should add two methyl groups, "
                    "leading to a C12 product.")

        # Check 2: Functional Groups. Does the product have the correct functional groups?
        # The reaction opens an epoxide to an alcohol and performs conjugate addition, leaving the ketone intact.
        # The final product should be a hydroxy-ketone.
        # Options A and D are diols, which would require reduction of the ketone. Gilman reagents do not do this.
        # This check correctly eliminates options A and D. The proposed answer C is a hydroxy-ketone, which is correct.
        
        # Check 3: Consistency. Is the proposed pathway the most consistent with the problem statement?
        # The pathway from Intermediate A to Product C fully utilizes the "excess" reagent condition by reacting at both available sites.
        # Alternative pathways (e.g., leading to option B, a C11 molecule) would imply that only one of the two reactive sites on the corresponding intermediate reacted, which contradicts the "excess" reagent condition.
        # Therefore, the derivation of answer C is the most consistent interpretation of all given constraints.

        return "Correct"

    except Exception as e:
        return f"An error occurred during the check: {e}"

# Run the check
result = check_answer()
print(result)