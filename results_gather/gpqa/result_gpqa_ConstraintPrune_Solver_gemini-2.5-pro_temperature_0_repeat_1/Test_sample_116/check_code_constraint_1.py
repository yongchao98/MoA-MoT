try:
    from rdkit import Chem
except ImportError:
    raise ImportError("This checking code requires the RDKit library. Please install it using 'pip install rdkit'")

def check_optical_isomerism(smiles_string: str) -> bool:
    """
    Checks if a molecule is chiral, which is the condition for optical isomerism.
    
    This function uses RDKit to parse a SMILES string and identify any potential
    stereocenters. This includes both classic tetrahedral chiral centers and other
    sources of chirality like atropisomerism (axial chirality), which is relevant
    for the biphenyl compound.
    
    Args:
        smiles_string: The SMILES representation of the molecule.
        
    Returns:
        True if the molecule is determined to be chiral, False otherwise.
    """
    mol = Chem.MolFromSmiles(smiles_string)
    if mol is None:
        # This would indicate an error in the SMILES string provided.
        return False
        
    # FindMolChiralCenters identifies potential stereocenters.
    # We set includeUnassigned=True to find potential centers even if their
    # specific R/S configuration isn't defined in the SMILES.
    # This function correctly identifies tetrahedral centers and atropisomers.
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    
    # If the list of chiral centers is not empty, the molecule is chiral.
    return len(chiral_centers) > 0

def check_answer():
    """
    Checks the correctness of the LLM's answer by analyzing each compound.
    """
    # Mapping of compound number to its name and SMILES string.
    # SMILES are a standard way to represent chemical structures.
    compounds = {
        1: {
            "name": "dimethyl 6,6'-dinitro-[1,1'-biphenyl]-2,2'-dicarboxylate",
            "smiles": "COC(=O)c1c(cccc1[N+](=O)[O-])c2c(C(=O)OC)cccc2[N+](=O)[O-]"
        },
        2: {
            "name": "methyl 2-hydroxypropanoate",
            "smiles": "CC(O)C(=O)OC"
        },
        3: {
            "name": "benzophenone",
            "smiles": "O=C(c1ccccc1)c1ccccc1"
        },
        4: {
            "name": "dimethyl fumarate",
            "smiles": "COC(=O)/C=C/C(=O)OC"  # The slashes define trans-geometry
        }
    }

    # The LLM's answer is B, which corresponds to compounds 1 and 2.
    llm_answer_indices = {1, 2}

    # Determine which compounds are chiral according to our programmatic check.
    code_determined_chiral_indices = set()
    for index, data in compounds.items():
        if check_optical_isomerism(data["smiles"]):
            code_determined_chiral_indices.add(index)
            
    # Compare the set of chiral compounds found by the code with the LLM's answer.
    if code_determined_chiral_indices == llm_answer_indices:
        return "Correct"
    else:
        # If they don't match, provide a detailed reason.
        reason = "The answer is incorrect.\n"
        reason += f"The code determined that the compounds showing optical isomerism are: {sorted(list(code_determined_chiral_indices))}.\n"
        reason += f"The provided answer claims the correct compounds are: {sorted(list(llm_answer_indices))}.\n"
        
        # Provide a breakdown of the analysis for clarity.
        reason += "\nAnalysis breakdown:\n"
        for i in range(1, 5):
            is_chiral = (i in code_determined_chiral_indices)
            is_in_answer = (i in llm_answer_indices)
            status = "Matches" if is_chiral == is_in_answer else "MISMATCH"
            reason += f" - Compound {i}: Code found chiral={is_chiral}, Answer implies chiral={is_in_answer}. ({status})\n"
            
        return reason

# Run the check and print the result.
result = check_answer()
print(result)