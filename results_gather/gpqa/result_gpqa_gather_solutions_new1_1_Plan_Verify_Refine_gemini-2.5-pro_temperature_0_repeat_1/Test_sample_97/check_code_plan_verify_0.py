try:
    from rdkit import Chem
    from rdkit.Chem import rdMolDescriptors
except ImportError:
    # If rdkit is not installed, we cannot run the check.
    # We will provide a message to the user.
    pass

def check_chemistry_answer():
    """
    Checks the correctness of the answer to the chemistry question by verifying
    the proposed reaction pathway using chemical principles.
    """
    try:
        Chem
    except NameError:
        return ("Could not perform check: The 'rdkit' library is required. "
                "Please install it using 'pip install rdkit-pypi'.")

    # --- Data Definition ---
    # The question maps letters to IUPAC names.
    # The final answer from the LLM is <<<A>>>.
    # We verify that option A, bicyclo[3.2.0]hept-6-ene, is the correct precursor.
    
    options = {
        "A": "bicyclo[3.2.0]hept-6-ene",
        "B": "2-methyl-3-methylenebicyclo[2.1.0]pentane",
        "C": "2-methylbicyclo[3.1.0]hex-2-ene",
        "D": "1,2-dimethylenecyclopentane"
    }
    
    # SMILES strings for all relevant molecules
    smiles_db = {
        "bicyclo[3.2.0]hept-6-ene": "C1=CC2C(C1)CCC2",
        "2-methyl-3-methylenebicyclo[2.1.0]pentane": "C=C1C(C)C2C1C2",
        "2-methylbicyclo[3.1.0]hex-2-ene": "CC1=CC2CC1C2",
        "1,2-dimethylenecyclopentane": "C=C1C(=C)CCC1",
        "1-propene": "C=CC",
        "1-(prop-1-en-1-yl)-2-vinylcyclopentane": "CCC=CC1C(C=C)CCC1"
    }

    llm_answer_key = "A"
    candidate_name = options[llm_answer_key]
    candidate_mol = Chem.MolFromSmiles(smiles_db[candidate_name])
    product_mol = Chem.MolFromSmiles(smiles_db["1-(prop-1-en-1-yl)-2-vinylcyclopentane"])

    # --- Verification Checks for the Chosen Answer ---

    # Check 1: Atom Conservation
    # The net reaction is A + 1-propene -> Product.
    # A (C7H10) + C3H6 -> C10H16
    formula_A = rdMolDescriptors.CalcMolFormula(candidate_mol)
    formula_product = rdMolDescriptors.CalcMolFormula(product_mol)
    if not (formula_A == "C7H10" and formula_product == "C10H16"):
        return f"Incorrect atom count. Candidate A is {formula_A} (expected C7H10) and product is {formula_product} (expected C10H16)."

    # Check 2: Must be a bicyclic alkene suitable for ROCM.
    if rdMolDescriptors.GetNumRings(candidate_mol) != 2:
        return f"Candidate '{candidate_name}' is not bicyclic. It has {rdMolDescriptors.GetNumRings(candidate_mol)} ring(s) and cannot undergo ROCM."

    # Check 3: Must preserve a cyclopentane core.
    # This means the double bond must be in the non-cyclopentane ring.
    sssr = Chem.GetSymmSSSR(candidate_mol)
    ring_sizes = [len(r) for r in sssr]
    if 5 not in ring_sizes:
        return f"Candidate '{candidate_name}' does not contain a 5-membered ring, so it cannot form a cyclopentane core."
    
    opened_ring_size = -1
    for bond in candidate_mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
            for ring in sssr:
                if bond.GetBeginAtomIdx() in ring and bond.GetEndAtomIdx() in ring:
                    opened_ring_size = len(ring)
                    break
            break
    
    if opened_ring_size == 5:
        return f"The double bond in '{candidate_name}' is in the 5-membered ring. ROCM would destroy the required cyclopentane core."

    # Check 4: Must produce a 1,2-disubstituted product.
    # This requires the bridgehead atoms to be adjacent.
    ring_info = candidate_mol.GetRingInfo()
    bridgehead_atoms = [i for i in range(candidate_mol.GetNumAtoms()) if ring_info.NumAtomRings(i) > 1]
    
    if len(bridgehead_atoms) != 2:
        return f"Candidate '{candidate_name}' is not a standard fused bicyclic system (found {len(bridgehead_atoms)} bridgehead atoms)."
        
    atom1_idx, atom2_idx = bridgehead_atoms[0], bridgehead_atoms[1]
    if not candidate_mol.GetBondBetweenAtoms(atom1_idx, atom2_idx):
        return f"The bridgehead atoms in '{candidate_name}' are not adjacent. This would lead to a 1,3- or 1,4-disubstituted product, not the required 1,2-disubstituted product."

    # --- Verify Rejection of Other Options ---
    # B) No 5-membered ring
    mol_b = Chem.MolFromSmiles(smiles_db[options["B"]])
    if 5 in [len(r) for r in Chem.GetSymmSSSR(mol_b)]:
        return "Logic check failed: Option B was expected to NOT have a 5-membered ring."
    
    # C) Double bond in 5-membered ring
    mol_c = Chem.MolFromSmiles(smiles_db[options["C"]])
    opened_ring_size_c = -1
    for bond in mol_c.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
            for ring in Chem.GetSymmSSSR(mol_c):
                if bond.GetBeginAtomIdx() in ring and bond.GetEndAtomIdx() in ring:
                    opened_ring_size_c = len(ring)
                    break
            break
    if opened_ring_size_c != 5:
        return "Logic check failed: Option C was expected to have its double bond in the 5-membered ring."

    # D) Not bicyclic
    mol_d = Chem.MolFromSmiles(smiles_db[options["D"]])
    if rdMolDescriptors.GetNumRings(mol_d) == 2:
        return "Logic check failed: Option D was expected to NOT be bicyclic."

    return "Correct"

# Run the check
result = check_chemistry_answer()
print(result)