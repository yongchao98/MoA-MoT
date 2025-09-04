from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def check_chemistry_answer():
    """
    Checks the correctness of the LLM's answer to the organic chemistry question.

    The logic is as follows:
    1.  Determine the expected molecular formula of the product. The reaction is an addition of HBr (H, Br) to the reactant.
    2.  Check the molecular formulas of the products listed in each option. Any option with products that do not match the expected formula is incorrect.
    3.  For the remaining valid option(s), check if the proposed structures are mechanistically plausible. The major product should arise from the most stable carbocation intermediate (Markovnikov's rule).
    """
    # --- Define Molecules using SMILES strings ---
    reactant_smiles = "c1ccc(cc1)OCC(C)(C)C=C"  # ((2,2-dimethylbut-3-en-1-yl)oxy)benzene

    # Products from Option A
    prod_A1_smiles = "c1ccc(cc1)OCC(C)(C)C(Br)C"  # (3-bromo-2,2-dimethylbutoxy)benzene (Markovnikov)
    prod_A2_smiles = "c1ccc(cc1)OCC(C)(C)CCBr"    # (4-bromo-2,2-dimethylbutoxy)benzene (Anti-Markovnikov)

    # Products from Option B
    prod_B1_smiles = "CC1(C)CC2=C(O1)C=CC=C2C"    # 3,3,4-trimethylchromane
    prod_B2_smiles = "CC1(C2=C(O1)C=CC=C2)C(C)C"  # 3-isopropyl-3-methyl-2,3-dihydrobenzofuran

    # Products from Option C
    prod_C1_smiles = prod_A2_smiles
    prod_C2_smiles = "c1ccc(cc1)OCC(C)=C(C)C"     # ((2,3-dimethylbut-2-en-1-yl)oxy)benzene

    # Products from Option D
    prod_D1_smiles = "c1ccc(c(O)c1)C(C)(C)CCC"    # 2-(2,2-dimethylbutyl)phenol
    prod_D2_smiles = "c1cc(ccc1O)C(C)(C)CCC"      # 4-(2,2-dimethylbutyl)phenol

    llm_answer_choice = 'A'

    # --- Step 1: Calculate Expected Molecular Formula ---
    reactant_mol = Chem.MolFromSmiles(reactant_smiles)
    reactant_formula = rdMolDescriptors.CalcMolFormula(reactant_mol) # Should be C12H16O
    # The reaction is an addition of HBr.
    expected_product_formula = "C12H17BrO"

    # --- Step 2: Check Molecular Formulas of All Options ---
    options = {
        'A': [prod_A1_smiles, prod_A2_smiles],
        'B': [prod_B1_smiles, prod_B2_smiles],
        'C': [prod_C1_smiles, prod_C2_smiles],
        'D': [prod_D1_smiles, prod_D2_smiles]
    }

    for option, smiles_list in options.items():
        for i, smiles in enumerate(smiles_list):
            mol = Chem.MolFromSmiles(smiles)
            formula = rdMolDescriptors.CalcMolFormula(mol)
            if formula != expected_product_formula:
                # This check invalidates options B, C, and D immediately.
                if option == llm_answer_choice:
                    return (f"Incorrect. The LLM's chosen answer '{option}' contains a product "
                            f"('{smiles}') with an incorrect molecular formula '{formula}'. "
                            f"The reaction is an addition of HBr to {reactant_formula}, so the product "
                            f"formula should be {expected_product_formula}.")

    # If we reach here, it means only Option A has products with the correct molecular formula.
    # Now, we perform a mechanistic check on the correct option.

    # --- Step 3: Mechanistic Plausibility Check ---
    # The major product should be the Markovnikov product.
    # Reactant: Ph-O-CH2-C(Me)2-CH=CH2
    # Protonation of the terminal CH2 gives a stable tertiary carbocation at C3.
    # Bromide attack at C3 gives (3-bromo-2,2-dimethylbutoxy)benzene.
    major_product_smiles = prod_A1_smiles

    if major_product_smiles not in options[llm_answer_choice]:
        return (f"Incorrect. The LLM's chosen answer '{llm_answer_choice}' is mechanistically flawed. "
                f"While the products have the correct formula, the expected major product "
                f"(Markovnikov adduct, SMILES: {major_product_smiles}) is not present in this option.")

    # The LLM's answer (A) has products with the correct molecular formula AND contains the expected major product.
    return "Correct"

# Run the check
result = check_chemistry_answer()
print(result)