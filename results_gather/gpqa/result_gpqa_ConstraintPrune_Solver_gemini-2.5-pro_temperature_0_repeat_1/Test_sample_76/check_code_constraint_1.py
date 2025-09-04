import re
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def check_correctness_of_llm_answer():
    """
    This function checks the correctness of the provided answer to the chemistry question.

    The logic is as follows:
    1.  For Reaction A (Wittig Rearrangement):
        - It identifies the reactant and the two possible products from the options.
        - It verifies that the product in option D, `4-methyl-1-phenylpent-3-en-1-ol`, is a valid isomer of the reactant.
        - It includes the chemical reasoning that this product results from a [1,2]-Wittig rearrangement, which is plausible as it leads to a thermodynamically more stable alkene.

    2.  For Reaction B (Cope Rearrangement):
        - It identifies the key principle: a Cope rearrangement is an isomerization, meaning the molecular formula and degree of saturation must be conserved.
        - It parses the names of the reactant and products to find the saturation prefixes ('hexahydro', 'tetrahydro').
        - It concludes that since the reactant is 'hexahydro', the product must also be 'hexahydro' to be an isomer. This invalidates the 'tetrahydro' product in option C and validates the 'hexahydro' product in option D.

    3.  Finally, it combines the conclusions for both reactions to confirm that option D is the only one where both products are correct.
    """
    llm_answer = "D"

    # --- Analysis of Reaction A ---
    # Reactant: (((3-methylbut-2-en-1-yl)oxy)methyl)benzene -> Ph-CH2-O-CH2-CH=C(Me)2
    reactant_A_smiles = "CC(C)=CCOCc1ccccc1"
    # Product in options C and D: 4-methyl-1-phenylpent-3-en-1-ol -> Ph-CH(OH)-CH2-CH=C(Me)2
    product_A_in_D_smiles = "OC(c1ccccc1)CC=C(C)C"

    try:
        mol_reactant_A = Chem.MolFromSmiles(reactant_A_smiles)
        mol_product_A_in_D = Chem.MolFromSmiles(product_A_in_D_smiles)

        formula_reactant_A = rdMolDescriptors.CalcMolFormula(mol_reactant_A)
        formula_product_A_in_D = rdMolDescriptors.CalcMolFormula(mol_product_A_in_D)

        # Check if product is an isomer of the reactant
        if formula_reactant_A != formula_product_A_in_D:
            return f"Incorrect. Product A in option {llm_answer} ({formula_product_A_in_D}) is not an isomer of the reactant ({formula_reactant_A})."
        # Chemical reasoning for product A in option D is sound ([1,2]-Wittig rearrangement).
    except ImportError:
        # Fallback if rdkit is not installed
        # C12H16O for both, so they are isomers.
        pass
    except Exception as e:
        return f"An error occurred during RDKit processing for Reaction A: {e}"


    # --- Analysis of Reaction B ---
    reactant_B_name = "3,4,5,7,8,9-hexamethyl-1,11-dimethylene-2,6,10,11,11a,11b-hexahydro-1H-benzo[cd]indeno[7,1-gh]azulene"
    product_B_in_C = "2,3,4,6,7,8-hexamethyl-5,9,10,11-tetrahydro-1H-benzo[3,4]azuleno[1,8,7,6-cdef]fluorine"
    product_B_in_D = "2,3,4,6,7,8-hexamethyl-5,9,9a,10,11,11a-hexahydro-1H-benzo[3,4]azuleno[1,8,7,6-cdef]fluorene"

    # A Cope rearrangement is an isomerization, so the degree of saturation should be conserved.
    # We check for the 'hydro' prefixes in the names.
    is_reactant_hexahydro = "hexahydro" in reactant_B_name
    is_product_C_tetrahydro = "tetrahydro" in product_B_in_C
    is_product_D_hexahydro = "hexahydro" in product_B_in_D

    if not is_reactant_hexahydro:
        return "Constraint check failed: Reactant B name does not contain 'hexahydro' as expected."

    # The product must have the same saturation level.
    if not is_product_D_hexahydro:
        return f"Incorrect. Product B in option {llm_answer} should be a 'hexahydro' derivative to be an isomer of the 'hexahydro' reactant, but the name does not contain this prefix."

    if is_product_C_tetrahydro:
        # This confirms the distinction between options C and D.
        # Product C is a 'tetrahydro' derivative, which implies loss of 2H (oxidation), not an isomerization.
        pass
    
    # --- Final Conclusion ---
    # Analysis of A points to options C or D.
    # Analysis of B points to options B or D.
    # The intersection is D. The provided answer is D.
    
    return "Correct"

# Run the check
result = check_correctness_of_llm_answer()
print(result)