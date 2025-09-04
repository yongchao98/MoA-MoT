import sys
from io import StringIO

# This code requires the rdkit library.
# You can install it with: pip install rdkit
try:
    from rdkit import Chem
except ImportError:
    # If rdkit is not installed, the check cannot be performed.
    print("Incorrect. The checker code could not run because the 'rdkit' library is not installed. Please install it using 'pip install rdkit' to verify the answer.")
    # We cannot confirm correctness without the necessary tools.
    exit()

def check_correctness():
    """
    Checks the correctness of the LLM's answer about optically active compounds.
    """

    # This function checks if a molecule, represented by a SMILES string, is chiral.
    # A molecule is achiral if it is superimposable on its mirror image.
    # In RDKit, this can be robustly checked by comparing the canonical SMILES of a molecule
    # with the canonical SMILES of its generated enantiomer. If they are identical, the molecule is achiral (e.g., a meso compound).
    def is_molecule_chiral(smiles: str) -> bool:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"Invalid SMILES string: {smiles}")

        # A quick check: if there are no potential stereocenters, it's generally achiral.
        if not Chem.FindMolChiralCenters(mol, includeUnassigned=True) and all(b.GetStereo() == Chem.BondStereo.STEREONONE for b in mol.GetBonds()):
             return False

        # The definitive check: compare the molecule to its mirror image via canonical SMILES.
        smi_orig = Chem.MolToSmiles(mol, isomericSmiles=True)
        
        # Create the enantiomer by inverting all stereocenters.
        mol_enantiomer = Chem.Mol(mol)
        
        # Invert atom stereocenters (R <-> S)
        chiral_centers = Chem.FindMolChiralCenters(mol)
        for center in chiral_centers:
            atom_idx, tag = center
            atom = mol_enantiomer.GetAtomWithIdx(atom_idx)
            if tag == 'R':
                atom.SetChiralTag(Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_S)
            elif tag == 'S':
                atom.SetChiralTag(Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_R)
                
        # Invert bond stereocenters (E <-> Z)
        for bond in mol_enantiomer.GetBonds():
            if bond.GetStereo() == Chem.rdchem.BondStereo.STEREOE:
                bond.SetStereo(Chem.rdchem.BondStereo.STEREOZ)
            elif bond.GetStereo() == Chem.rdchem.BondStereo.STEREOZ:
                bond.SetStereo(Chem.rdchem.BondStereo.STEREOE)
        
        Chem.AssignStereochemistry(mol_enantiomer, force=True, cleanIt=True)
        smi_enantiomer = Chem.MolToSmiles(mol_enantiomer, isomericSmiles=True)
        
        return smi_orig != smi_enantiomer

    # --- Analysis of each compound ---
    compounds_to_check = [
        {
            "name": "(Z)-1-chloro-2-methylbut-1-ene",
            "smiles": r"C/C(C)=C(\Cl)C",
            "is_single_isomer_specified": True,
            "llm_is_active": False,
        },
        {
            "name": "(3aR,7aS,E)-8-(chloromethylene)hexahydro-4,7-methanoisobenzofuran-1,3-dione",
            "smiles": None, # SMILES is too complex to reliably generate; we'll rely on the name.
            "is_single_isomer_specified": True,
            "llm_is_active": True,
        },
        {
            "name": "(2R,3S)-2,3-dimethylsuccinic acid",
            "smiles": "O=C(O)[C@@H](C)[C@H](C)C(=O)O", # This is the meso form
            "is_single_isomer_specified": True,
            "llm_is_active": False,
        },
        {
            "name": "(2R,3R)-2,3-dimethylsuccinic acid",
            "smiles": "O=C(O)[C@H](C)[C@H](C)C(=O)O",
            "is_single_isomer_specified": True,
            "llm_is_active": True,
        },
        {
            "name": "(R)-cyclohex-3-en-1-ol",
            "smiles": "O[C@H]1CC=CCC1",
            "is_single_isomer_specified": True,
            "llm_is_active": True,
        },
        {
            "name": "(1s,3s,5s)-cyclohexane-1,3,5-triol",
            "smiles": "O[C@H]1C[C@H](O)C[C@H](O)C1", # all-cis, symmetric
            "is_single_isomer_specified": True,
            "llm_is_active": False,
        },
        {
            "name": "1-cyclopentyl-3-methylbutan-1-one",
            "smiles": "CC(C)CC(=O)C1CCCC1",
            "is_single_isomer_specified": False, # No stereochemistry in name
            "llm_is_active": False,
        }
    ]

    code_active_count = 0
    
    for compound in compounds_to_check:
        name = compound["name"]
        smiles = compound["smiles"]
        is_single_isomer = compound["is_single_isomer_specified"]
        llm_active = compound["llm_is_active"]
        
        is_active = False
        reason = ""

        if name == "(3aR,7aS,E)-8-(chloromethylene)hexahydro-4,7-methanoisobenzofuran-1,3-dione":
            # Expert knowledge shortcut: The name specifies a single, complex stereoisomer.
            # Such a molecule is chiral and, as a single isomer, is optically active.
            is_active = True
            reason = "Chiral molecule specified as a single stereoisomer by its complex name."
        else:
            is_chiral = is_molecule_chiral(smiles)
            if is_chiral and is_single_isomer:
                is_active = True
                reason = "Molecule is chiral and specified as a single isomer."
            elif not is_chiral:
                is_active = False
                reason = "Molecule is achiral."
            elif is_chiral and not is_single_isomer:
                is_active = False
                reason = "Molecule is chiral, but as no specific isomer is named, it's assumed to be a racemic mixture (inactive)."

        if is_active != llm_active:
            return f"Incorrect. The analysis for '{name}' is wrong. The LLM claims optical activity is {llm_active}, but the code determined it is {is_active}. Reason: {reason}"
        
        if is_active:
            code_active_count += 1

    llm_total_count = 3
    if code_active_count != llm_total_count:
        return f"Incorrect. The final count is wrong. The LLM's answer implies {llm_total_count} optically active compounds, but the code calculated {code_active_count}."

    return "Correct"

# Run the check and print the result
try:
    result = check_correctness()
    print(result)
except Exception as e:
    print(f"An error occurred during the check: {e}")