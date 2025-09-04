# First, ensure you have RDKit installed:
# pip install rdkit

import sys

def check_correctness():
    """
    This function checks the correctness of the given answer by analyzing the chirality
    of each compound using the RDKit library.
    """
    try:
        from rdkit import Chem
    except ImportError:
        print("This checking code requires the RDKit library. Please install it using 'pip install rdkit'")
        return

    # The answer provided by the LLM, converted to an integer.
    # D corresponds to 3.
    llm_answer_map = {'A': 4, 'B': 5, 'C': 2, 'D': 3}
    llm_answer_key = "D"
    llm_answer_count = llm_answer_map.get(llm_answer_key)

    if llm_answer_count is None:
        print(f"Invalid answer key '{llm_answer_key}'. Must be one of A, B, C, D.")
        return

    # A list of dictionaries, each containing the compound's name and its
    # SMILES representation or a special analysis rule.
    compounds_to_check = [
        {
            "name": "(Z)-1-chloro-2-methylbut-1-ene",
            "smiles": "Cl/C=C(\\C)CC", # Z-isomer, but no chiral center
            "analysis": "rdkit"
        },
        {
            "name": "(3aR,7aS,E)-8-(chloromethylene)hexahydro-4,7-methanoisobenzofuran-1,3-dione",
            "smiles": None, # Structure is too complex for simple SMILES generation
            "analysis": "by_name" # The name specifies a single enantiomer (3aR,7aS)
        },
        {
            "name": "(2R,3S)-2,3-dimethylsuccinic acid",
            "smiles": "C[C@H]([C@@H](C)C(=O)O)C(=O)O", # This is the meso form
            "analysis": "rdkit"
        },
        {
            "name": "(2R,3R)-2,3-dimethylsuccinic acid",
            "smiles": "C[C@H]([C@H](C)C(=O)O)C(=O)O", # This is a chiral enantiomer
            "analysis": "rdkit"
        },
        {
            "name": "(R)-cyclohex-3-en-1-ol",
            "smiles": "O[C@H]1CC=CCC1", # Chiral center at C1
            "analysis": "rdkit"
        },
        {
            "name": "(1s,3s,5s)-cyclohexane-1,3,5-triol",
            "smiles": "O[C@H]1C[C@H](O)C[C@H](O)C1", # all-cis isomer, which is achiral due to symmetry
            "analysis": "rdkit"
        },
        {
            "name": "1-cyclopentyl-3-methylbutan-1-one",
            "smiles": "CC(C)CC(=O)C1CCCC1", # No chiral centers
            "analysis": "rdkit"
        }
    ]

    active_compounds_list = []
    analysis_log = []

    for compound in compounds_to_check:
        name = compound["name"]
        is_active = False
        reason = ""

        if compound["analysis"] == "by_name":
            # Heuristic: If a specific enantiomer is named with R/S descriptors, it is chiral.
            is_active = True
            reason = "Considered optically active because the name specifies a single enantiomer (e.g., R/S)."
        
        elif compound["analysis"] == "rdkit":
            mol = Chem.MolFromSmiles(compound["smiles"])
            # Generate canonical SMILES with stereochemistry.
            # A molecule is optically active if it's chiral. A robust proxy for this in RDKit
            # is the presence of '@' in the canonical isomeric SMILES. This correctly identifies
            # meso compounds as achiral because their '@' symbols are removed during canonicalization.
            canon_smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
            if '@' in canon_smiles:
                is_active = True
                reason = f"Optically active. Its canonical SMILES ('{canon_smiles}') retains chiral center markers."
            else:
                is_active = False
                reason = f"Not optically active. Its canonical SMILES ('{canon_smiles}') has no chiral center markers, indicating it's achiral."

        if is_active:
            active_compounds_list.append(name)
        
        analysis_log.append(f"  - {name}: {'Active' if is_active else 'Inactive'}\n    Reason: {reason}")

    calculated_count = len(active_compounds_list)

    if calculated_count == llm_answer_count:
        print("Correct")
    else:
        error_message = (
            f"Incorrect. The provided answer is {llm_answer_count}, but the analysis indicates there are {calculated_count} optically active compounds.\n\n"
            f"Detailed Analysis:\n"
        )
        error_message += "\n".join(analysis_log)
        error_message += f"\n\nList of optically active compounds found ({calculated_count}):\n"
        if active_compounds_list:
            for c in active_compounds_list:
                error_message += f"- {c}\n"
        else:
            error_message += "- None\n"
            
        print(error_message)

check_correctness()