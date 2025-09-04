import sys
try:
    from rdkit import Chem
except ImportError:
    print("RDKit is not installed. Please install it using 'pip install rdkit'")
    # Exit gracefully if rdkit is not available
    sys.exit("RDKit dependency not met.")

def is_chiral(smiles: str) -> bool:
    """
    Checks if a molecule represented by a SMILES string is chiral.
    This function correctly handles achiral molecules (including meso compounds)
    by comparing the canonical SMILES of a molecule with its enantiomer.

    Args:
        smiles: The SMILES string of the molecule.

    Returns:
        True if the molecule is chiral, False otherwise.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        # Invalid SMILES string cannot be processed.
        raise ValueError(f"Invalid SMILES: {smiles}")

    # If a molecule has no potential chiral centers, it is achiral.
    # This is a quick check for simple achiral molecules.
    if not Chem.FindMolChiralCenters(mol, includeUnassigned=True):
        return False

    # The robust check: A molecule is chiral if and only if it is not
    # superimposable on its mirror image. In RDKit, this means the
    # canonical SMILES of a molecule and its enantiomer will be different.
    
    # Canonicalize the original molecule's SMILES
    canon_smi_original = Chem.MolToSmiles(mol, isomericSmiles=True)

    # If RDKit's canonicalization removes all stereochemistry, it's achiral (meso).
    if '@' not in canon_smi_original:
        return False

    # Create the enantiomer by inverting all stereocenters.
    # Note: This creates a new molecule instance.
    enantiomer = Chem.MolFromSmiles(canon_smi_original)
    Chem.InvertMol(enantiomer)
    
    # Canonicalize the enantiomer's SMILES
    canon_smi_enantiomer = Chem.MolToSmiles(enantiomer, isomericSmiles=True)

    # If the canonical representations are different, the molecule is chiral.
    return canon_smi_original != canon_smi_enantiomer

def check_answer():
    """
    Checks the correctness of the LLM's answer by analyzing each compound.
    """
    # Dictionary mapping compound names to their SMILES representation.
    # A special value 'CHIRAL_BY_NAME' is used for compounds whose chirality
    # is guaranteed by their specific IUPAC name.
    compounds_data = {
        "(Z)-1-chloro-2-methylbut-1-ene": "Cl/C=C(\\C)CC",
        "(3aR,7aS,E)-8-(chloromethylene)hexahydro-4,7-methanoisobenzofuran-1,3-dione": "CHIRAL_BY_NAME",
        "(2R,3S)-2,3-dimethylsuccinic acid": "C[C@H]([C@@H](C)C(=O)O)C(=O)O",
        "(2R,3R)-2,3-dimethylsuccinic acid": "C[C@H]([C@H](C)C(=O)O)C(=O)O",
        "(R)-cyclohex-3-en-1-ol": "O[C@H]1CC=CCC1",
        "(1s,3s,5s)-cyclohexane-1,3,5-triol": "O[C@H]1C[C@H](O)C[C@H](O)C1",
        "1-cyclopentyl-3-methylbutan-1-one": "CC(C)CC(=O)C1CCCC1"
    }

    # The list of compounds the provided answer claims are optically active.
    llm_active_compounds = {
        "(3aR,7aS,E)-8-(chloromethylene)hexahydro-4,7-methanoisobenzofuran-1,3-dione",
        "(2R,3R)-2,3-dimethylsuccinic acid",
        "(R)-cyclohex-3-en-1-ol"
    }
    
    # The total count claimed by the answer.
    llm_count = 3

    code_active_compounds = set()
    errors = []

    for name, data in compounds_data.items():
        is_active = False
        if data == "CHIRAL_BY_NAME":
            # The name contains specific stereodescriptors (R/S), indicating a single
            # enantiomer, which is by definition optically active.
            is_active = True
        else:
            try:
                is_active = is_chiral(data)
            except Exception as e:
                errors.append(f"Error processing '{name}': {e}")
                continue
        
        if is_active:
            code_active_compounds.add(name)

    # Check 1: Compare the total count of active compounds.
    if len(code_active_compounds) != llm_count:
        return (f"Incorrect. The answer claims there are {llm_count} optically active compounds, "
                f"but the code found {len(code_active_compounds)}. "
                f"The compounds determined to be active are: {list(code_active_compounds)}")

    # Check 2: Compare the set of active compounds.
    if code_active_compounds != llm_active_compounds:
        missing_from_llm = code_active_compounds - llm_active_compounds
        extra_in_llm = llm_active_compounds - code_active_compounds
        error_msg = "Incorrect. The set of optically active compounds does not match.\n"
        if missing_from_llm:
            error_msg += f"The answer failed to identify the following as active: {list(missing_from_llm)}\n"
        if extra_in_llm:
            error_msg += f"The answer incorrectly identified the following as active: {list(extra_in_llm)}\n"
        return error_msg

    if errors:
        return f"Execution completed with errors: {errors}"

    return "Correct"

# Run the check and print the result.
result = check_answer()
print(result)