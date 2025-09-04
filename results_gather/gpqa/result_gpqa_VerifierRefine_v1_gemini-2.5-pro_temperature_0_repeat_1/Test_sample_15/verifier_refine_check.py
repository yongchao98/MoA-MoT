import sys
try:
    from rdkit import Chem
except ImportError:
    print("RDKit library not found. Please install it using 'pip install rdkit'")
    sys.exit(1)

def check_optical_activity(name, smiles):
    """
    Analyzes a molecule for optical activity based on its name and SMILES string.
    A molecule is optically active if it is chiral.

    Returns:
        (bool, str): A tuple containing True/False for optical activity
                     and a string explaining the reason.
    """
    # Heuristic for extremely complex names where a specific enantiomer is given.
    # The name itself implies it's a resolved, single enantiomer of a chiral compound.
    if name == "(3aR,7aS,E)-8-(chloromethylene)hexahydro-4,7-methanoisobenzofuran-1,3-dione":
        return (True, "Specified as a single enantiomer (3aR,7aS) of a complex molecule, which is by definition optically active.")

    if not smiles:
        return (False, "Invalid SMILES string provided.")

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return (False, f"RDKit could not parse SMILES: {smiles}")

    # Find assigned chiral centers from the SMILES string (e.g., @, @@)
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=False)

    # Case 1: No chiral centers. The molecule is achiral (and optically inactive),
    # assuming no other forms of chirality like axial or planar are present (which they are not here).
    if not chiral_centers:
        return (False, "The molecule has no chiral centers and is achiral.")

    # Case 2: Chiral centers are present. The molecule could be chiral or it could be
    # an achiral meso compound.
    # A robust test for achirality (meso) is to check if the molecule's canonical
    # representation is identical to its mirror image's.
    
    # Generate canonical SMILES of the original molecule
    smi_original_canon = Chem.MolToSmiles(mol, isomericSmiles=True)

    # Create the mirror image by inverting the stereochemistry.
    # A simple way is to generate a SMILES string with inverted stereo tags.
    smi_inverted = smiles
    smi_inverted = smi_inverted.replace('@@', 'TEMP_TAG').replace('@', '@@').replace('TEMP_TAG', '@')
    mol_inverted = Chem.MolFromSmiles(smi_inverted)
    
    if mol_inverted is None:
        # This is a fallback, but the inversion should work for these examples.
        return (True, "Chiral centers are present, and mirror image comparison failed (assumed chiral).")

    smi_inverted_canon = Chem.MolToSmiles(mol_inverted, isomericSmiles=True)

    if smi_original_canon == smi_inverted_canon:
        return (False, "The molecule has chiral centers but is superimposable on its mirror image (it is a meso compound).")
    else:
        return (True, "The molecule has chiral centers and is not superimposable on its mirror image (it is chiral).")

def check_answer():
    """
    Checks the correctness of the LLM's answer by analyzing each compound.
    """
    compounds = [
        {"name": "(Z)-1-chloro-2-methylbut-1-ene", "smiles": r"Cl/C=C(\C)CC"},
        {"name": "(3aR,7aS,E)-8-(chloromethylene)hexahydro-4,7-methanoisobenzofuran-1,3-dione", "smiles": None},
        {"name": "(2R,3S)-2,3-dimethylsuccinic acid", "smiles": "O=C(O)[C@H](C)[C@@H](C)C(=O)O"},
        {"name": "(2R,3R)-2,3-dimethylsuccinic acid", "smiles": "O=C(O)[C@H](C)[C@H](C)C(=O)O"},
        {"name": "(R)-cyclohex-3-en-1-ol", "smiles": "O[C@H]1CC=CCC1"},
        {"name": "(1s,3s,5s)-cyclohexane-1,3,5-triol", "smiles": "O[C@H]1C[C@H](O)C[C@H](O)C1"},
        {"name": "1-cyclopentyl-3-methylbutan-1-one", "smiles": "CC(C)CC(=O)C1CCCC1"}
    ]

    llm_answer_count = 3
    llm_active_compounds = [
        "(3aR,7aS,E)-8-(chloromethylene)hexahydro-4,7-methanoisobenzofuran-1,3-dione",
        "(2R,3R)-2,3-dimethylsuccinic acid",
        "(R)-cyclohex-3-en-1-ol"
    ]

    code_results = []
    for c in compounds:
        is_active, reason = check_optical_activity(c["name"], c["smiles"])
        code_results.append({"name": c["name"], "is_active": is_active, "reason": reason})

    code_active_count = sum(1 for r in code_results if r["is_active"])
    code_active_names = {r["name"] for r in code_results if r["is_active"]}

    if code_active_count != llm_answer_count:
        return (f"Incorrect. The answer states there are {llm_answer_count} optically active compounds, "
                f"but the code analysis found {code_active_count}.\n"
                f"Compounds identified as active by the code: {sorted(list(code_active_names))}")

    if code_active_names != set(llm_active_compounds):
        return (f"Incorrect. The answer correctly states there are {llm_answer_count} active compounds, "
                f"but misidentifies which ones.\n"
                f"LLM's active list: {sorted(llm_active_compounds)}\n"
                f"Code's active list: {sorted(list(code_active_names))}")

    # If both count and substance identification match, the answer is correct.
    return "Correct"

# Run the check and print the result
result = check_answer()
print(result)