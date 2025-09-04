import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors

def check_cycloaddition_answer():
    """
    Checks the correctness of the answer for the Diels-Alder reaction.
    The verification process involves three main steps:
    1.  Check the bridge atom type based on the reactants.
    2.  Check the stereochemistry (EXO vs. ENDO) based on the reaction conditions (Heat).
    3.  Verify that the provided answer corresponds to the correct structure using computational chemistry.
    """
    # --- Problem Definition ---
    # Question: Identify the EXO product of the [4+2] cycloaddition reaction.
    # Reactants: 2,5-dimethylthiophene + Furan-2,5-dione + Heat
    # Provided Answer from the other LLM: B
    
    options = {
        "A": "(3aR,4S,7R,7aS)-4,7-dimethyl-3a,4,7,7a-tetrahydro-4,7-epithioisobenzofuran-1,3-dione",
        "B": "(3aR,4R,7S,7aS)-4,7-dimethyl-3a,4,7,7a-tetrahydro-4,7-epithioisobenzofuran-1,3-dione",
        "C": "(3aR,4R,7S,7aS)-4,7-dimethyl-3a,4,7,7a-tetrahydro-4,7-epoxybenzo[c]thiophene-1,3-dione",
        "D": "(3aR,4S,7R,7aS)-4,7-dimethyl-3a,4,7,7a-tetrahydro-4,7-epoxybenzo[c]thiophene-1,3-dione"
    }
    llm_answer = "B"

    # --- Constraint 1: Verify the Bridge Atom ---
    # The diene is 2,5-dimethylthiophene. The sulfur atom from the thiophene ring
    # becomes the bridge atom in the bicyclic product.
    # Therefore, the product must be an "epithio" adduct, not an "epoxy" (oxygen) adduct.
    
    if "epoxy" in options[llm_answer]:
        return f"Incorrect. The answer {llm_answer} is wrong because the product must be an 'epithio' adduct, not an 'epoxy' adduct. The sulfur atom from the thiophene diene forms the bridge."
    
    # This check confirms that C and D are incorrect, which aligns with the provided reasoning.

    # --- Constraint 2: Verify Stereochemistry (EXO vs. ENDO) ---
    # The question specifies "Heat", which favors the thermodynamically more stable product.
    # In Diels-Alder reactions, the EXO adduct is generally more stable than the ENDO adduct due to less steric hindrance.
    # Therefore, we need to identify which option (A or B) represents the EXO adduct.
    
    # The SMILES strings from the provided answer are used for verification.
    # These correspond to the IUPAC names for the ENDO (A) and EXO (B) diastereomers.
    smiles_map = {
        "A": "C[C@@]12C=C[C@](C)(S1)[C@H]1C(=O)O[C@H]1C2=O",  # Claimed ENDO
        "B": "C[C@@]12C=C[C@](C)(S1)[C@@H]1C(=O)O[C@@H]1C2=O" # Claimed EXO
    }

    # Sanity check: Verify molecular formula for the given answer's SMILES
    try:
        mol_b = Chem.MolFromSmiles(smiles_map[llm_answer])
        if mol_b is None:
            return f"Error: The SMILES string for the given answer '{llm_answer}' is invalid."
        formula = Descriptors.rdMolDescriptors.CalcMolFormula(mol_b)
        expected_formula = "C10H10O3S" # C6H8S + C4H2O3 = C10H10O3S
        if formula != expected_formula:
            return f"Incorrect. The molecular formula for the product should be {expected_formula}, but the SMILES for option {llm_answer} corresponds to {formula}."
    except KeyError:
        return f"Error: No SMILES string provided for the answer option {llm_answer}."

    # Perform geometric analysis to identify the true EXO isomer.
    identified_exo_option = None
    for option, smiles in smiles_map.items():
        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(mol)
        
        # Generate a 3D conformer
        params = AllChem.ETKDGv3()
        params.randomSeed = 42  # for reproducibility
        if AllChem.EmbedMolecule(mol, params) == -1:
            return f"Error: RDKit failed to generate a 3D conformer for option {option}, cannot verify stereochemistry."
        
        conf = mol.GetConformer()
        
        try:
            # Identify key atoms for the geometric check
            s_idx = mol.GetSubstructMatch(Chem.MolFromSmarts('[S]'))[0]
            s_atom = mol.GetAtomWithIdx(s_idx)
            bridgehead_indices = [a.GetIdx() for a in s_atom.GetNeighbors() if a.GetAtomicNum() == 6]
            
            # Find the two carbons of the anhydride ring that are part of the fusion
            anhydride_match = mol.GetSubstructMatch(Chem.MolFromSmarts('C1C(=O)OC(=O)C1'))
            fusion_indices = [idx for idx in anhydride_match if mol.GetAtomWithIdx(idx).GetDegree() > 2]

            if len(bridgehead_indices) != 2 or len(fusion_indices) != 2:
                raise ValueError("Failed to identify key atoms for geometric analysis.")

            # Calculate positions and vectors
            s_pos = np.array(conf.GetAtomPosition(s_idx))
            # A centroid of the 4 atoms at the base of the bicyclic system
            base_centroid = (np.array(conf.GetAtomPosition(bridgehead_indices[0])) + 
                             np.array(conf.GetAtomPosition(bridgehead_indices[1])) + 
                             np.array(conf.GetAtomPosition(fusion_indices[0])) + 
                             np.array(conf.GetAtomPosition(fusion_indices[1]))) / 4.0
            
            # A point representing the anhydride part of the molecule
            anhydride_midpoint = (np.array(conf.GetAtomPosition(fusion_indices[0])) + np.array(conf.GetAtomPosition(fusion_indices[1]))) / 2.0

            # Vector from the base to the sulfur bridge
            vec_to_s = s_pos - base_centroid
            # Vector from the base to the anhydride part
            vec_to_anhydride = anhydride_midpoint - base_centroid
            
            # The dot product determines the relative orientation.
            # Negative dot product means the vectors point in opposite directions (trans), which is EXO.
            # Positive dot product means they point in the same direction (cis), which is ENDO.
            dot_product = np.dot(vec_to_s, vec_to_anhydride)
            
            if dot_product < 0:
                identified_exo_option = option

        except Exception as e:
            return f"Geometric check failed for option {option}: {e}"

    # --- Constraint 3: Final Verification ---
    # The question asks for the EXO product.
    if identified_exo_option is None:
        return "Incorrect. The geometric check could not identify an EXO product among options A and B. The provided SMILES may both be ENDO or there is an error in the analysis."
        
    if identified_exo_option == llm_answer:
        return "Correct"
    else:
        return f"Incorrect. The question asks for the EXO product, which is favored by heat. The geometric analysis shows that option '{identified_exo_option}' is the EXO adduct, but the provided answer was '{llm_answer}'."

# Execute the check and print the result
result = check_cycloaddition_answer()
print(result)