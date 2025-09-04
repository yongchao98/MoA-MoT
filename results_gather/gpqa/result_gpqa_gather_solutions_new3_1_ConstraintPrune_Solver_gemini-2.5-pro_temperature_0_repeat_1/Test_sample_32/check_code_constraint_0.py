import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np
import warnings

# Suppress RDKit warnings for cleaner output
from rdkit import rdBase
rdBase.DisableLog('rdApp.*')
warnings.filterwarnings("ignore", category=UserWarning, module='pubchempy')

def check_cycloaddition_product():
    """
    Checks the correctness of the answer by verifying nomenclature and stereochemistry.
    """
    # Step 1: Define the options and the provided answer to check
    options = {
        'A': "(3aR,4R,7S,7aS)-4,7-dimethyl-3a,4,7,7a-tetrahydro-4,7-epoxybenzo[c]thiophene-1,3-dione",
        'B': "(3aR,4R,7S,7aS)-4,7-dimethyl-3a,4,7,7a-tetrahydro-4,7-epithioisobenzofuran-1,3-dione",
        'C': "(3aR,4S,7R,7aS)-4,7-dimethyl-3a,4,7,7a-tetrahydro-4,7-epoxybenzo[c]thiophene-1,3-dione",
        'D': "(3aR,4S,7R,7aS)-4,7-dimethyl-3a,4,7,7a-tetrahydro-4,7-epithioisobenzofuran-1,3-dione",
    }
    provided_answer_key = 'D'

    # Step 2: Constraint 1 - Check nomenclature for "epithio" bridge and "isobenzofuran" core
    print("Constraint 1: Checking product nomenclature...")
    valid_candidates = {}
    eliminated_candidates = {}
    for key, name in options.items():
        if "epithio" in name and "isobenzofuran" in name:
            valid_candidates[key] = name
        else:
            reason = "Incorrect nomenclature: should be 'epithioisobenzofuran'."
            if "epoxy" in name:
                reason = "Incorrect bridge: should be 'epithio' (sulfur), not 'epoxy' (oxygen)."
            eliminated_candidates[key] = reason
    
    print(f"  - Eliminated: {list(eliminated_candidates.keys())} ({eliminated_candidates})")
    print(f"  - Valid candidates based on name: {list(valid_candidates.keys())}")

    if provided_answer_key not in valid_candidates:
        return f"Incorrect: The provided answer '{provided_answer_key}' is eliminated by the nomenclature check. {eliminated_candidates.get(provided_answer_key)}"

    # Step 3: Constraint 2 & 3 - Identify the EXO isomer via geometric check
    print("\nConstraint 2 & 3: Identifying the EXO isomer...")
    print("  - Principle: 'Heat' favors the thermodynamically stable EXO product.")
    print("  - Method: Checking geometry of candidates B and D.")

    exo_candidate = None
    endo_candidate = None

    for key, name in valid_candidates.items():
        try:
            # Use PubChem to get a SMILES string from the IUPAC name
            compounds = pcp.get_compounds(name, 'name')
            if not compounds:
                print(f"  - Warning: Could not find candidate '{key}' in PubChem. Cannot perform geometric check.")
                continue
            
            mol_smiles = compounds[0].isomeric_smiles
            mol = Chem.MolFromSmiles(mol_smiles)
            mol = Chem.AddHs(mol)

            # Generate a 3D conformer
            if AllChem.EmbedMolecule(mol, randomSeed=42) == -1:
                print(f"  - Warning: RDKit failed to generate a 3D conformer for candidate '{key}'.")
                continue
            
            conf = mol.GetConformer()

            # Find key atoms to define geometry
            s_idx = mol.GetSubstructMatch(Chem.MolFromSmarts('[s]'))[0]
            anhydride_match = mol.GetSubstructMatch(Chem.MolFromSmarts('C1C(=O)OC(=O)C1'))
            bridgehead_C_indices = [a.GetIdx() for a in mol.GetAtomWithIdx(s_idx).GetNeighbors() if a.GetAtomicNum() == 6]
            fusion_C_indices = [idx for idx in anhydride_match if mol.GetAtomWithIdx(idx).GetDegree() > 2]

            # Calculate vectors from the ring centroid to the S-bridge and the anhydride ring
            all_bridge_carbons_pos = [np.array(conf.GetAtomPosition(i)) for i in bridgehead_C_indices + fusion_C_indices]
            centroid = np.mean(all_bridge_carbons_pos, axis=0)
            s_pos = np.array(conf.GetAtomPosition(s_idx))
            fusion_midpoint = np.mean([np.array(conf.GetAtomPosition(i)) for i in fusion_C_indices], axis=0)
            
            vec_to_s = s_pos - centroid
            vec_to_anhydride = fusion_midpoint - centroid
            dot_product = np.dot(vec_to_s, vec_to_anhydride)

            # EXO: anhydride is anti to S-bridge (vectors point opposite, dot product < 0)
            # ENDO: anhydride is syn to S-bridge (vectors point same way, dot product > 0)
            if dot_product < 0:
                print(f"  - Candidate '{key}' has EXO geometry (anti).")
                exo_candidate = key
            else:
                print(f"  - Candidate '{key}' has ENDO geometry (syn).")
                endo_candidate = key

        except Exception as e:
            return f"An error occurred during the geometric check for candidate {key}: {e}"

    # Step 4: Final verification
    print("\n--- Conclusion ---")
    if exo_candidate is None:
        return "Incorrect: The geometric check failed to identify an EXO product among the valid candidates."
        
    if exo_candidate == provided_answer_key:
        return "Correct"
    else:
        return f"Incorrect: The provided answer is '{provided_answer_key}', but the geometric analysis identified candidate '{exo_candidate}' as the EXO product. The provided answer corresponds to the ENDO product."

# Run the check
result = check_cycloaddition_product()
print(f"\nFinal Result: {result}")