import sys

def check_answer():
    """
    Checks the correctness of the LLM's answer by analyzing the chemical properties
    of the given compounds using the RDKit library.
    """
    try:
        from rdkit import Chem
    except ImportError:
        return "Error: RDKit library not found. Please install it using 'pip install rdkit-pypi' to run this check."

    # --- Part A: Tautomerism Analysis ---
    # Keto-enol tautomerism requires an alpha-hydrogen on an sp3 carbon.
    # We will check for this structural feature.

    def has_keto_enol_tautomerism_site(smiles: str) -> bool:
        """
        Checks for the presence of an sp3 alpha-hydrogen next to a carbonyl group.
        This is a strong indicator of keto-enol tautomerism.
        """
        mol = Chem.MolFromSmiles(smiles)
        if not mol: return False
        
        # Pattern for a carbonyl group with an adjacent sp3 carbon having at least one hydrogen.
        # [C;H1,H2,H3][CX3](=O) means: a Carbon with 1, 2, or 3 Hydrogens,
        # next to a trigonal Carbon double-bonded to an Oxygen.
        tautomer_pattern = Chem.MolFromSmarts('[C;H1,H2,H3][CX3](=O)')
        return mol.HasSubstructMatch(tautomer_pattern)

    # Define compounds for Part A
    compounds_A = {
        "benzoquinone": "O=C1C=CC(=O)C=C1",
        "cyclohexane-1,3,5-trione": "O=C1CC(=O)CC(=O)C1"
    }

    # Find the compound that does NOT show tautomerism
    compound_A_result = None
    for name, smiles in compounds_A.items():
        if not has_keto_enol_tautomerism_site(smiles):
            compound_A_result = name
            
    # --- Part B: Optical Isomerism Analysis ---
    # Optical isomerism requires chirality. We check for chiral centers.

    def has_optical_isomerism(smiles: str) -> bool:
        """
        Checks if a molecule has chiral centers, a common cause of optical isomerism.
        """
        mol = Chem.MolFromSmiles(smiles)
        if not mol: return False
        
        # Find chiral centers. includeUnassigned=True is important for structures
        # where chirality isn't explicitly marked in the SMILES string.
        chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
        return len(chiral_centers) > 0

    # Define compounds for Part B
    compounds_B = {
        "methyl 2-hydroxypropanoate": "CC(O)C(=O)OC",
        "dimethyl fumarate": "COC(=O)/C=C/C(=O)OC"
    }

    # Find the compound that DOES show optical isomerism
    compound_B_result = None
    for name, smiles in compounds_B.items():
        if has_optical_isomerism(smiles):
            compound_B_result = name

    # --- Verification ---
    # The provided answer corresponds to option B.
    llm_answer_A = "benzoquinone"
    llm_answer_B = "methyl 2-hydroxypropanoate"

    errors = []
    if compound_A_result != llm_answer_A:
        reason = (f"The compound that does not show tautomerism is '{compound_A_result}', not '{llm_answer_A}'. "
                  f"Benzoquinone lacks sp3 alpha-hydrogens, while cyclohexane-1,3,5-trione has them, allowing it to tautomerize to the stable aromatic phloroglucinol.")
        errors.append(f"Incorrect analysis for Part A. {reason}")

    if compound_B_result != llm_answer_B:
        reason = (f"The compound that shows optical isomerism is '{compound_B_result}', not '{llm_answer_B}'. "
                  f"Methyl 2-hydroxypropanoate has a chiral center (a carbon bonded to -H, -OH, -CH3, and -COOCH3), while dimethyl fumarate is an achiral planar molecule.")
        errors.append(f"Incorrect analysis for Part B. {reason}")

    if not errors:
        return "Correct"
    else:
        return "\n".join(errors)

# Run the check and print the result
result = check_answer()
print(result)