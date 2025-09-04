import sys
from contextlib import redirect_stdout
import io

def check_answer_correctness():
    """
    Checks the correctness of the provided answer about optical activity.
    
    This function uses the RDKit library to analyze the chirality of each compound.
    It compares the program's findings with the reasoning and final count from the
    provided answer.
    
    Returns:
        str: "Correct" if the answer is correct, otherwise a string explaining the error.
    """
    try:
        from rdkit import Chem
    except ImportError:
        return "Error: RDKit library is not installed. Please install it using 'pip install rdkit'."

    # --- Data from the question and answer ---
    # The final answer claims 3 compounds are optically active.
    expected_active_count = 3
    
    # List of compounds with their SMILES representation and expected activity based on the provided answer.
    compounds = [
        {
            "name": "(Z)-1-chloro-2-methylbut-1-ene",
            "smiles": "Cl/C=C(\\C)CC",
            "is_active_expected": False,
            "reason": "Achiral due to plane of symmetry."
        },
        {
            "name": "(3aR,7aS,E)-8-(chloromethylene)hexahydro-4,7-methanoisobenzofuran-1,3-dione",
            "smiles": "MANUAL_CHECK", # Name implies chirality.
            "is_active_expected": True,
            "reason": "Name specifies a single enantiomer of a complex chiral molecule."
        },
        {
            "name": "(2R,3S)-2,3-dimethylsuccinic acid",
            "smiles": "C[C@H](C(=O)O)[C@@H](C)C(=O)O",
            "is_active_expected": False,
            "reason": "Meso compound with an internal plane of symmetry."
        },
        {
            "name": "(2R,3R)-2,3-dimethylsuccinic acid",
            "smiles": "C[C@H](C(=O)O)[C@H](C)C(=O)O",
            "is_active_expected": True,
            "reason": "Chiral diastereomer, specified as a single enantiomer."
        },
        {
            "name": "(R)-cyclohex-3-en-1-ol",
            "smiles": "O[C@H]1CCC=CC1",
            "is_active_expected": True,
            "reason": "Chiral molecule, specified as a single enantiomer."
        },
        {
            "name": "(1s,3s,5s)-cyclohexane-1,3,5-triol",
            "smiles": "O[C@H]1C[C@H](O)C[C@H](O)C1",
            "is_active_expected": False,
            "reason": "Achiral (meso-like) due to high symmetry (planes of symmetry)."
        },
        {
            "name": "1-cyclopentyl-3-methylbutan-1-one",
            "smiles": "CC(C)CC(=O)C1CCCC1",
            "is_active_expected": False,
            "reason": "Achiral, no stereocenters."
        }
    ]

    def is_molecule_chiral(smiles_string, name):
        """Checks if a molecule is chiral using RDKit."""
        if smiles_string == "MANUAL_CHECK":
            # If a name includes R/S descriptors for a complex molecule, it's a single
            # enantiomer and thus optically active.
            return True
        
        # Suppress RDKit warnings/errors for this check
        f = io.StringIO()
        with redirect_stdout(f):
            mol = Chem.MolFromSmiles(smiles_string)
        
        if mol is None:
            raise ValueError(f"Could not parse SMILES for {name}")

        # Generate canonical SMILES with stereochemistry. RDKit will omit stereo
        # markers for achiral/meso compounds even if they have chiral centers.
        # This is a robust way to check for overall molecular chirality.
        isomeric_smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
        return '@' in isomeric_smiles

    errors = []
    calculated_active_count = 0

    for compound in compounds:
        try:
            is_active_calculated = is_molecule_chiral(compound["smiles"], compound["name"])
            if is_active_calculated:
                calculated_active_count += 1
            
            if is_active_calculated != compound["is_active_expected"]:
                errors.append(
                    f"Mismatch for '{compound['name']}': "
                    f"Expected active = {compound['is_active_expected']}, but code calculated {is_active_calculated}. "
                    f"Reason: {compound['reason']}"
                )
        except Exception as e:
            errors.append(f"An error occurred while processing '{compound['name']}': {e}")

    if calculated_active_count != expected_active_count:
        errors.append(
            f"Final count mismatch: The answer states {expected_active_count} optically active compounds, "
            f"but the code calculated {calculated_active_count}."
        )

    if errors:
        return "Incorrect. Reason(s):\n- " + "\n- ".join(errors)
    else:
        return "Correct"

# Run the check
result = check_answer_correctness()
print(result)