import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

def check_diels_alder_product():
    """
    This function checks the correctness of the provided LLM answer for a Diels-Alder reaction.

    The check involves three main steps:
    1.  **Chemical Plausibility**: Verifies that the product type matches the reactants.
        The reaction of a thiophene (containing sulfur) should yield an "epithio" bridge,
        not an "epoxy" bridge (containing oxygen).
    2.  **Thermodynamic Plausibility**: Verifies the stereochemical outcome based on reaction conditions.
        The condition "Heat" drives the reaction towards the more stable thermodynamic product,
        which in a Diels-Alder reaction is the EXO adduct.
    3.  **Geometric Verification**: Checks if the structure given for the proposed answer
        actually has the required EXO geometry. This is done by generating a 3D model and
        analyzing the relative positions of the sulfur bridge and the anhydride ring.
    """
    # --- Step 1 & 2: Check Chemical and Thermodynamic Plausibility ---
    llm_answer_choice = "A"
    options = {
        "A": "(3aR,4S,7R,7aS)-4,7-dimethyl-3a,4,7,7a-tetrahydro-4,7-epithioisobenzofuran-1,3-dione",
        "B": "(3aR,4R,7S,7aS)-4,7-dimethyl-3a,4,7,7a-tetrahydro-4,7-epithioisobenzofuran-1,3-dione",
        "C": "(3aR,4R,7S,7aS)-4,7-dimethyl-3a,4,7,7a-tetrahydro-4,7-epoxybenzo[c]thiophene-1,3-dione",
        "D": "(3aR,4S,7R,7aS)-4,7-dimethyl-3a,4,7,7a-tetrahydro-4,7-epoxybenzo[c]thiophene-1,3-dione"
    }

    # Check if the answer correctly identifies an "epithio" product.
    if "epoxy" in options[llm_answer_choice]:
        return "Incorrect. The selected answer is an 'epoxy' derivative. The reaction involves 2,5-dimethylthiophene, which contains a sulfur atom, so the product must be an 'epithio' derivative with a sulfur bridge."

    # The question asks for the EXO product, which is correct for a reaction under "Heat".
    # The LLM claims that option A is the EXO product. We must verify this specific claim.

    # --- Step 3: Geometric Verification of the LLM's chosen structure ---
    # The LLM's response provided a SMILES string for its answer 'A'. We will analyze it.
    llm_smiles_for_A = "C[C@@]12C=C[C@](C)(S1)[C@H]1C(=O)O[C@H](C1=O)C2"

    mol = Chem.MolFromSmiles(llm_smiles_for_A)
    if not mol:
        return f"Incorrect. The SMILES string '{llm_smiles_for_A}' provided for answer A in the LLM's explanation is invalid and cannot be parsed by RDKit."

    mol = Chem.AddHs(mol)
    if AllChem.EmbedMolecule(mol, randomSeed=42) == -1:
        # If embedding fails, the structure is likely highly strained or invalid.
        return "Error in checker: RDKit failed to generate a valid 3D conformer for the provided SMILES of answer A. This suggests the SMILES string may represent an impossible structure."

    try:
        AllChem.UFFOptimizeMolecule(mol)
    except Exception:
        # A failure to optimize is not critical, the initial embedding is usually sufficient.
        pass

    conf = mol.GetConformer()

    try:
        # Identify key atoms for defining the molecule's geometry
        s_idx = mol.GetSubstructMatch(Chem.MolFromSmarts('[S]'))[0]
        bridgehead_indices = [a.GetIdx() for a in mol.GetAtomWithIdx(s_idx).GetNeighbors() if a.GetAtomicNum() == 6]

        anhydride_pattern = Chem.MolFromSmarts('C1C(=O)OC(=O)C1')
        anhydride_match = mol.GetSubstructMatch(anhydride_pattern)
        fusion_indices = [idx for idx in anhydride_match if mol.GetAtomWithIdx(idx).GetDegree() > 2 and mol.GetAtomWithIdx(idx).GetAtomicNum() == 6]

        anhydride_o_pattern = Chem.MolFromSmarts('C(=O)OC(=O)')
        anhydride_o_idx = mol.GetSubstructMatch(anhydride_o_pattern)[1]

        if len(bridgehead_indices) != 2 or len(fusion_indices) != 2:
            raise ValueError("Failed to identify key structural atoms (2 bridgeheads and 2 fusion carbons).")

    except Exception as e:
        return f"Error in checker: Could not identify the necessary atoms for geometric analysis. This may be due to an unexpected molecular structure from the SMILES string. Reason: {e}"

    # Calculate vectors from a central point to the sulfur bridge and the anhydride ring
    s_pos = np.array(conf.GetAtomPosition(s_idx))
    anhydride_o_pos = np.array(conf.GetAtomPosition(anhydride_o_idx))

    # Define a reference point as the centroid of the four base carbons of the bicyclic system
    base_centroid = (np.array(conf.GetAtomPosition(bridgehead_indices[0])) +
                     np.array(conf.GetAtomPosition(bridgehead_indices[1])) +
                     np.array(conf.GetAtomPosition(fusion_indices[0])) +
                     np.array(conf.GetAtomPosition(fusion_indices[1]))) / 4

    # Vector from the centroid to the sulfur bridge
    vec_to_s = s_pos - base_centroid
    # Vector from the centroid to the anhydride's ether oxygen (points "out" towards the anhydride ring)
    vec_to_anhydride_o = anhydride_o_pos - base_centroid

    # The dot product of these vectors determines the relative orientation.
    # If negative, vectors point in opposite directions -> S and anhydride are TRANS -> EXO
    # If positive, vectors point in the same direction -> S and anhydride are CIS -> ENDO
    dot_product = np.dot(vec_to_s, vec_to_anhydride_o)

    is_exo = dot_product < 0

    # Final Verdict based on the geometric analysis
    if is_exo:
        # The LLM's claim that structure A is EXO is consistent with our analysis.
        # Since the EXO product is expected, the answer is correct.
        return "Correct"
    else:
        # The LLM's claim is contradicted by our analysis.
        return f"Incorrect. The question asks for the EXO product, which is favored by heat. The LLM claims option A is the EXO product. However, a geometric analysis of the structure provided for option A shows it has an ENDO configuration, where the sulfur bridge and the anhydride ring are on the same side (cis) of the molecule's base (dot product = {dot_product:.2f})."

# Run the check
result = check_diels_alder_product()
print(result)