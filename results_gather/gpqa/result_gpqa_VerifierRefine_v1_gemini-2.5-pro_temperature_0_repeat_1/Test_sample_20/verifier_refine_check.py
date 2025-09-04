# First, ensure you have RDKit installed:
# pip install rdkit-pypi

try:
    from rdkit import Chem
except ImportError:
    # This block allows the code to run in environments where rdkit is not pre-installed.
    # It provides a clear error message if the necessary library is missing.
    print("RDKit library not found. Please install it to run this check: pip install rdkit-pypi")
    # We will exit if rdkit is not available, as the check cannot proceed.
    # In a real application, you might return a specific error code or message.
    # For this purpose, we'll raise the error to stop execution.
    raise

def check_keto_enol_tautomerism(mol: Chem.Mol) -> bool:
    """
    Checks for the structural requirement for keto-enol tautomerism.
    It looks for an sp3-hybridized alpha-carbon with at least one hydrogen.
    An alpha-carbon is a carbon atom adjacent to a carbonyl carbon.
    """
    # Pattern for a carbonyl group (C=O)
    carbonyl_pattern = Chem.MolFromSmarts('[CX3]=[OX1]')
    if not mol.HasSubstructMatch(carbonyl_pattern):
        return False

    # Iterate through all atoms to find alpha-carbons with hydrogens
    for atom in mol.GetAtoms():
        # We are looking for sp3 carbons with at least one hydrogen
        if atom.GetAtomicNum() == 6 and \
           atom.GetHybridization() == Chem.rdchem.HybridizationType.SP3 and \
           atom.GetTotalNumHs() > 0:
            # Check if this carbon is an alpha-carbon (neighbor to a carbonyl carbon)
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 6 and \
                   neighbor.GetHybridization() == Chem.rdchem.HybridizationType.SP2 and \
                   any(b.GetBondType() == Chem.rdchem.BondType.DOUBLE and b.GetOtherAtom(neighbor).GetAtomicNum() == 8 for b in neighbor.GetBonds()):
                    # This is an sp3 alpha-carbon with a hydrogen. Tautomerism is possible.
                    return True
    return False

def check_optical_isomerism(mol: Chem.Mol) -> bool:
    """
    Checks for optical isomerism by finding potential chiral centers.
    A molecule with at least one chiral center is typically optically active.
    """
    # FindMolChiralCenters identifies atoms that are chiral centers.
    # includeUnassigned=True ensures it finds potential centers even if stereochemistry isn't specified in the input SMILES.
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    return len(chiral_centers) > 0

def check_answer():
    """
    Main function to check the correctness of the LLM's answer.
    """
    # Define molecules using their SMILES strings
    compounds = {
        "benzoquinone": Chem.MolFromSmiles("O=C1C=CC(=O)C=C1"),
        "cyclohexane-1,3,5-trione": Chem.MolFromSmiles("O=C1CC(=O)CC(=O)C1"),
        "methyl 2-hydroxypropanoate": Chem.MolFromSmiles("CC(O)C(=O)OC"),
        "dimethyl fumarate": Chem.MolFromSmiles("C/C(=C/C(=O)OC)/C(=O)OC")
    }

    # --- Part A: Tautomerism ---
    # The question asks for the compound that DOES NOT show tautomerism.
    shows_tautomerism = {
        "benzoquinone": check_keto_enol_tautomerism(compounds["benzoquinone"]),
        "cyclohexane-1,3,5-trione": check_keto_enol_tautomerism(compounds["cyclohexane-1,3,5-trione"])
    }

    # Based on our check, the compound that does NOT show tautomerism is the one with a 'False' result.
    compound_A = None
    for name, can_tautomerize in shows_tautomerism.items():
        if not can_tautomerize:
            compound_A = name
            break
    
    # --- Part B: Optical Isomerism ---
    # The question asks for the compound that WILL show optical isomerism.
    shows_optical_isomerism = {
        "methyl 2-hydroxypropanoate": check_optical_isomerism(compounds["methyl 2-hydroxypropanoate"]),
        "dimethyl fumarate": check_optical_isomerism(compounds["dimethyl fumarate"])
    }

    # Based on our check, the compound that WILL show optical isomerism is the one with a 'True' result.
    compound_B = None
    for name, is_chiral in shows_optical_isomerism.items():
        if is_chiral:
            compound_B = name
            break

    # The LLM's answer is D, which corresponds to:
    # A = benzoquinone
    # B = methyl 2-hydroxypropanoate
    llm_answer_A = "benzoquinone"
    llm_answer_B = "methyl 2-hydroxypropanoate"

    # --- Verification ---
    if compound_A != llm_answer_A:
        return f"Incorrect. For Part A, the compound that does not show tautomerism was determined to be '{compound_A}', but the answer states it is '{llm_answer_A}'. Benzoquinone lacks sp3 alpha-hydrogens, so it cannot tautomerize. Cyclohexane-1,3,5-trione has them and can."
    
    if compound_B != llm_answer_B:
        return f"Incorrect. For Part B, the compound that shows optical isomerism was determined to be '{compound_B}', but the answer states it is '{llm_answer_B}'. Methyl 2-hydroxypropanoate has a chiral center, while dimethyl fumarate is achiral."

    return "Correct"

# Run the check and print the result
result = check_answer()
print(result)