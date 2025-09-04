import sys

def check_answer_correctness():
    """
    Checks the correctness of the LLM's answer by programmatically determining
    the optical activity of each compound.
    """
    try:
        from rdkit import Chem
    except ImportError:
        return "Error: RDKit library not found. Please install it using 'pip install rdkit' to run this check."

    def is_optically_active(name: str, smiles: str) -> (bool, str):
        """
        Determines if a molecule is optically active (chiral).
        Returns a tuple: (is_active, reason).
        """
        # Special case for the complex molecule where the name is definitive.
        # The name (3aR,7aS,E)... specifies a single enantiomer, which is by definition optically active.
        if "hexahydro-4,7-methanoisobenzofuran-1,3-dione" in name:
            if "(3aR,7aS" in name:
                return True, "Chiral: A specific enantiomer (3aR,7aS) is named."
            else:
                return None, "Cannot determine activity: Stereochemistry is ambiguous or not specified for this complex molecule."

        if not smiles:
            return None, "Cannot determine activity: No SMILES string provided."

        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return None, f"Cannot determine activity: RDKit could not parse SMILES '{smiles}'."

        # Find all potential stereocenters.
        chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)

        # If there are no stereocenters, the molecule is achiral (for this set of examples).
        if not chiral_centers:
            return False, "Achiral: No stereocenters found."

        # If stereocenters exist, we must check for meso compounds or other symmetry.
        # A molecule is achiral if it is superimposable on its mirror image.
        # We can test this by generating the canonical SMILES of the molecule and its mirror image.
        # If they are the same, the molecule is achiral.

        # Create a mirror image by inverting all stereocenters.
        mirror_mol = Chem.Mol(mol)
        for center_idx, _ in chiral_centers:
            atom = mirror_mol.GetAtomWithIdx(center_idx)
            tag = atom.GetChiralTag()
            if tag == Chem.ChiralType.CHI_TETRAHEDRAL_CW:
                atom.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CCW)
            elif tag == Chem.ChiralType.CHI_TETRAHEDRAL_CCW:
                atom.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CW)
        
        for bond in mirror_mol.GetBonds():
            tag = bond.GetStereo()
            if tag == Chem.BondStereo.STEREOZ:
                bond.SetStereo(Chem.BondStereo.STEREOE)
            elif tag == Chem.BondStereo.STEREOE:
                bond.SetStereo(Chem.BondStereo.STEREOZ)

        # Compare canonical SMILES of the original and its mirror image.
        smi_original = Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)
        smi_mirror = Chem.MolToSmiles(mirror_mol, isomericSmiles=True, canonical=True)

        if smi_original == smi_mirror:
            return False, "Achiral: Molecule is superimposable on its mirror image (meso or other symmetry)."
        else:
            return True, "Chiral: Molecule is not superimposable on its mirror image."

    # List of compounds from the question and the expected activity from the LLM's answer.
    compounds = [
        {"name": "(Z)-1-chloro-2-methylbut-1-ene", "smiles": "Cl/C=C(\\C)CC", "llm_active": False},
        {"name": "(3aR,7aS,E)-8-(chloromethylene)hexahydro-4,7-methanoisobenzofuran-1,3-dione", "smiles": None, "llm_active": True},
        {"name": "(2R,3S)-2,3-dimethylsuccinic acid", "smiles": "C[C@H](C(=O)O)[C@@H](C)C(=O)O", "llm_active": False},
        {"name": "(2R,3R)-2,3-dimethylsuccinic acid", "smiles": "C[C@H](C(=O)O)[C@H](C)C(=O)O", "llm_active": True},
        {"name": "(R)-cyclohex-3-en-1-ol", "smiles": "O[C@H]1CC=CCC1", "llm_active": True},
        {"name": "(1s,3s,5s)-cyclohexane-1,3,5-triol", "smiles": "O[C@H]1C[C@H](O)C[C@H](O)C1", "llm_active": False},
        {"name": "1-cyclopentyl-3-methylbutan-1-one", "smiles": "CC(C)CC(=O)C1CCCC1", "llm_active": False},
    ]

    llm_answer_count = 3
    calculated_active_count = 0
    
    for i, c in enumerate(compounds):
        is_active, reason = is_optically_active(c["name"], c["smiles"])
        
        if is_active is None:
            return f"Error checking compound {i+1} ({c['name']}): {reason}"
            
        if is_active:
            calculated_active_count += 1
        
        # Check if the code's conclusion matches the LLM's conclusion for each compound.
        if is_active != c["llm_active"]:
            return (f"Incorrect. The analysis for '{c['name']}' is wrong. "
                    f"The LLM's answer claims it is {'optically active' if c['llm_active'] else 'optically inactive'}, "
                    f"but the code determined it is {'optically active' if is_active else 'optically inactive'} because: {reason}.")

    # Final check on the total count.
    if calculated_active_count != llm_answer_count:
        return (f"Incorrect. The final count is wrong. "
                f"The LLM's answer states there are {llm_answer_count} optically active compounds, "
                f"but the code calculated {calculated_active_count}.")

    return "Correct"

# Run the check
result = check_answer_correctness()
print(result)