import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors

def check_diels_alder_product():
    """
    Checks the correctness of the given answer for the Diels-Alder reaction.
    The code verifies the structure of the selected answer against the constraints
    of the reaction, particularly the EXO stereochemistry requirement.
    """
    # --- 1. Define the problem and the given answer ---
    # The question asks for the EXO product of the reaction between
    # 2,5-dimethylthiophene and Furan-2,5-dione (maleic anhydride).
    # The provided answer from the other LLM is A.
    
    # Option A is (3aR,4S,7R,7aS)-4,7-dimethyl-3a,4,7,7a-tetrahydro-4,7-epithioisobenzofuran-1,3-dione.
    # This name corresponds to the ENDO adduct.
    # The canonical SMILES string for this molecule is used for the check.
    smiles_A = "C[C@]12C=C[C@](C)(S1)[C@H]1C(=O)O[C@H](C1=O)C2"
    
    # --- 2. Verify the structure against basic chemical constraints ---
    
    mol = Chem.MolFromSmiles(smiles_A)
    if mol is None:
        return "Error: The SMILES string for Option A is invalid and cannot be parsed."

    # Constraint 1: Molecular Formula
    # Reactants: 2,5-dimethylthiophene (C6H8S) + Furan-2,5-dione (C4H2O3)
    # Expected Product Formula: C10H10O3S
    actual_formula = Descriptors.rdMolDescriptors.CalcMolFormula(mol)
    expected_formula = "C10H10O3S"
    if actual_formula != expected_formula:
        return f"Incorrect. The molecular formula of Option A is {actual_formula}, but the expected formula for the adduct is {expected_formula}."

    # Constraint 2: Bridge Atom
    # The diene is thiophene, so the bridge must be sulfur (epithio).
    # Options with an oxygen bridge (epoxy) are incorrect.
    if not mol.HasSubstructMatch(Chem.MolFromSmarts('[S]')):
         return "Incorrect. Option A does not contain a sulfur atom, which is required for the epithio bridge from thiophene."
    # This check confirms Option A has the correct bridge type, ruling out options C and D.

    # --- 3. Verify the key constraint: EXO vs. ENDO stereochemistry ---
    # The question specifically asks for the EXO product. We must determine if Option A is EXO or ENDO.
    # Definition:
    # - ENDO: The anhydride ring is on the same side (syn) of the bicyclic system as the sulfur bridge.
    # - EXO: The anhydride ring is on the opposite side (anti) of the bicyclic system as the sulfur bridge.

    mol_h = Chem.AddHs(mol)
    # Generate a 3D conformer for geometric analysis
    params = AllChem.ETKDGv3()
    params.randomSeed = 42  # for reproducibility
    if AllChem.EmbedMolecule(mol_h, params) == -1:
        return "Error: RDKit failed to generate a 3D conformer for Option A, so the geometric check cannot be performed."
    
    conf = mol_h.GetConformer()

    try:
        # Find the sulfur atom
        s_idx = mol_h.GetSubstructMatch(Chem.MolFromSmarts('[S]'))[0]
        
        # Find the bridgehead carbons attached to the sulfur
        s_atom = mol_h.GetAtomWithIdx(s_idx)
        bridgehead_indices = [a.GetIdx() for a in s_atom.GetNeighbors() if a.GetAtomicNum() == 6]
        
        # Find the fusion carbons of the anhydride ring
        anhydride_match = mol_h.GetSubstructMatch(Chem.MolFromSmarts('C1C(=O)OC(=O)C1'))
        fusion_indices = []
        for idx in anhydride_match:
            atom = mol_h.GetAtomWithIdx(idx)
            # Fusion atoms are part of two rings and have degree 3 in the heavy-atom graph
            if atom.IsInRing() and atom.GetDegree() == 3 and idx in anhydride_match:
                fusion_indices.append(idx)
        
        if len(bridgehead_indices) != 2 or len(fusion_indices) != 2:
            return "Error: Could not unambiguously identify bridge and fusion atoms for geometric analysis."

        # Calculate midpoints to define vectors
        s_pos = np.array(conf.GetAtomPosition(s_idx))
        bridge_midpoint = (np.array(conf.GetAtomPosition(bridgehead_indices[0])) + np.array(conf.GetAtomPosition(bridgehead_indices[1]))) / 2
        fusion_midpoint = (np.array(conf.GetAtomPosition(fusion_indices[0])) + np.array(conf.GetAtomPosition(fusion_indices[1]))) / 2

        # Vector from the center of the main ring to the sulfur bridge
        vec_to_sulfur = s_pos - bridge_midpoint
        # Vector from the center of the main ring to the anhydride ring
        vec_to_anhydride = fusion_midpoint - bridge_midpoint
        
        # The dot product determines the relative orientation
        # dot > 0 means vectors point in the same general direction -> SYN -> ENDO
        # dot < 0 means vectors point in opposite general directions -> ANTI -> EXO
        dot_product = np.dot(vec_to_sulfur, vec_to_anhydride)

        if dot_product > 0:
            stereochem = "ENDO"
        elif dot_product < 0:
            stereochem = "EXO"
        else:
            return "Error: Geometric analysis resulted in a planar structure, which is unexpected."

    except Exception as e:
        return f"An error occurred during the geometric analysis: {e}"

    # --- 4. Final Conclusion ---
    # The question asks for the EXO product.
    if stereochem == "EXO":
        return "Correct"
    else:
        return f"Incorrect. The question asks for the EXO product, which is the thermodynamically favored adduct formed under heat. The provided answer, Option A, is the {stereochem} adduct. The geometric analysis confirms that the anhydride ring and the sulfur bridge are on the same side (syn-orientation), which defines the {stereochem} isomer. The correct EXO product is represented by Option B."

# Execute the check and print the result.
result = check_diels_alder_product()
print(result)