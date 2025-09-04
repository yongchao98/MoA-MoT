try:
    from rdkit import Chem
    from rdkit.Chem import rdMolDescriptors
except ImportError:
    raise ImportError("This checking code requires the RDKit library. Please install it using 'pip install rdkit-pypi'")

def check_correctness():
    """
    Checks the correctness of the LLM's answer by analyzing both reactions.
    """
    # --- Part 1: Analysis of Reaction A ([2,3]-Wittig Rearrangement) ---

    # Reactant: (((3-methylbut-2-en-1-yl)oxy)methyl)benzene -> Prenyl benzyl ether
    reactant_A_smiles = "CC(C)=CCOCc1ccccc1"

    # Product A from the selected option C: 4-methyl-1-phenylpent-3-en-1-ol
    llm_product_A_smiles = "CC(C)=CCC(O)c1ccccc1"

    # The mechanistically correct product of the [2,3]-Wittig rearrangement
    # is 2,2-dimethyl-1-phenylbut-3-en-1-ol.
    correct_product_A_smiles = "C=CC(C)(C)C(O)c1ccccc1"

    # Convert SMILES to canonical SMILES for a robust comparison
    # We ignore stereochemistry as it's not the primary point of error here.
    mol_llm_product = Chem.MolFromSmiles(llm_product_A_smiles)
    mol_correct_product = Chem.MolFromSmiles(correct_product_A_smiles)

    if not mol_llm_product or not mol_correct_product:
        return "Failed to parse SMILES strings for Reaction A products."

    canon_llm_product = Chem.MolToSmiles(mol_llm_product, isomericSmiles=False)
    canon_correct_product = Chem.MolToSmiles(mol_correct_product, isomericSmiles=False)

    # Check if the LLM's chosen product matches the correct product
    is_product_A_correct = (canon_llm_product == canon_correct_product)

    # --- Part 2: Analysis of Reaction B (Cope Rearrangement) ---

    # The LLM's logic is that a Cope rearrangement is an isomerization,
    # so the degree of saturation ("hexahydro") must be conserved.
    reactant_B_name = "3,4,5,7,8,9-hexamethyl-1,11-dimethylene-2,6,10,11,11a,11b-hexahydro-1H-benzo[cd]indeno[7,1-gh]azulene"
    product_B_name_in_C = "2,3,4,6,7,8-hexamethyl-5,9,9a,10,11,11a-hexahydro-1H-benzo[3,4]azuleno[1,8,7,6-cdef]fluorene"
    product_B_name_in_B = "2,3,4,6,7,8-hexamethyl-5,9,10,11-tetrahydro-1H-benzo[3,4]azuleno[1,8,7,6-cdef]fluorine"

    is_logic_B_correct = ("hexahydro" in reactant_B_name and
                          "hexahydro" in product_B_name_in_C and
                          "tetrahydro" in product_B_name_in_B)

    # --- Final Verdict ---
    if is_product_A_correct and is_logic_B_correct:
        return "Correct"
    
    if not is_product_A_correct:
        error_reason = (
            "The answer is incorrect because the reasoning for Reaction A is flawed.\n"
            f"The LLM identifies product A as '4-methyl-1-phenylpent-3-en-1-ol' (SMILES: {llm_product_A_smiles}).\n"
            "However, the established mechanism for the [2,3]-Wittig rearrangement of the starting material "
            f"yields '2,2-dimethyl-1-phenylbut-3-en-1-ol' (SMILES: {correct_product_A_smiles}).\n"
            "The product selected by the LLM has a different molecular structure from the actual product. "
            "While the LLM's logic for Reaction B (conservation of saturation) is sound, its incorrect conclusion for Reaction A makes its final answer invalid."
        )
        return error_reason
        
    # This case is unlikely but included for completeness
    if not is_logic_B_correct:
        return "The answer is incorrect because the reasoning for Reaction B is flawed."

# Execute the check and print the result
result = check_correctness()
print(result)
