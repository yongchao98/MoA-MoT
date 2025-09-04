import sys
import io
from contextlib import redirect_stdout

# Suppress RDKit warnings for cleaner output
try:
    from rdkit import rdBase
    rdBase.DisableLog('rdApp.warning')
except ImportError:
    pass

def check_chemistry_answer():
    """
    Checks the correctness of the LLM's answer by analyzing the stereochemical
    features of the provided SMILES strings against the expected reaction outcome.
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
        import numpy as np
    except ImportError:
        return "Error: RDKit library not found. Please install it with 'pip install rdkit'."

    # --- Data from the question ---
    options = {
        'A': "O=C(OC)[C@@H]1[C@@H](C=C2)[C@@H]3[C@@H]([C@H]4[C@@H]3[C@H]5C=C[C@@H]4C5)[C@@H]2[C@@H]1C(OC)=O",
        'B': "O=C(OC)[C@H]1[C@@H](C=C2)[C@@H]3[C@@H]([C@@H]4[C@H]3[C@@H]5C=C[C@H]4C5)[C@@H]2[C@H]1C(OC)=O",
        'C': "O=C(OC)[C@@H]1[C@@H](C=C2)[C@@H]3[C@@H]([C@@H]4[C@H]3[C@H]5C=C[C@@H]4C5)[C@@H]2[C@H]1C(OC)=O",
        'D': "O=C(OC)[C@@H]1[C@H](C=C2)[C@@H]3[C@@H]([C@H]4[C@@H]3[C@@H]5C=C[C@H]4C5)[C@H]2[C@@H]1C(OC)=O"
    }
    llm_answer = 'A'

    # --- Analysis Functions ---

    def check_ester_config(smiles: str) -> str:
        """Checks if the two ester groups are cis or trans using RDKit."""
        try:
            mol = Chem.MolFromSmiles(smiles)
            if not mol: return "Invalid SMILES"
            
            patt = Chem.MolFromSmarts('[C;H1](C(=O)OC)')
            matches = mol.GetSubstructMatches(patt)
            
            if len(matches) != 2:
                return f"Error: Found {len(matches)} ester groups, expected 2."
            
            c1_idx, c2_idx = matches[0][0], matches[1][0]
            bond = mol.GetBondBetweenAtoms(c1_idx, c2_idx)
            
            if not bond:
                return "Error: Ester-bearing carbons are not adjacent."
            
            stereo = bond.GetStereo()
            if stereo == Chem.BondStereo.STEREOCIS:
                return "cis"
            elif stereo == Chem.BondStereo.STEREOTRANS:
                return "trans"
            else:
                return "undetermined"
        except Exception as e:
            return f"Error checking ester config: {e}"

    def check_addition_config(smiles: str) -> str:
        """Checks if the cyclopentadiene addition was syn or anti to the esters."""
        try:
            mol = Chem.MolFromSmiles(smiles)
            if not mol: return "Invalid SMILES"

            patt_ester = Chem.MolFromSmarts('[C;H1](C(=O)OC)')
            ester_matches = mol.GetSubstructMatches(patt_ester)
            if len(ester_matches) != 2: return f"Error: Found {len(ester_matches)} ester groups."
            ester_indices = [m[0] for m in ester_matches]

            patt_norbornane = Chem.MolFromSmarts('C1CC2CCC1C2')
            norbornane_matches = mol.GetSubstructMatches(patt_norbornane)
            if not norbornane_matches: return "Error: Norbornane substructure not found."
            bridge_idx = norbornane_matches[0][-1]
            
            bridge_atom = mol.GetAtomWithIdx(bridge_idx)
            attach_indices = [n.GetIdx() for n in bridge_atom.GetNeighbors()]

            mol_h = Chem.AddHs(mol)
            if AllChem.EmbedMolecule(mol_h, randomSeed=42) == -1:
                return "Error: 3D embedding failed."
            conf = mol_h.GetConformer()

            pos_ester1 = np.array(conf.GetAtomPosition(ester_indices[0]))
            pos_ester2 = np.array(conf.GetAtomPosition(ester_indices[1]))
            pos_bridge = np.array(conf.GetAtomPosition(bridge_idx))
            pos_attach1 = np.array(conf.GetAtomPosition(attach_indices[0]))
            pos_attach2 = np.array(conf.GetAtomPosition(attach_indices[1]))

            centroid_attach = (pos_attach1 + pos_attach2) / 2.0
            v_bridge = pos_bridge - centroid_attach
            centroid_ester = (pos_ester1 + pos_ester2) / 2.0
            v_ester = centroid_ester - centroid_attach

            dot_product = np.dot(v_bridge, v_ester)
            
            if dot_product > 0.1: return "syn"
            elif dot_product < -0.1: return "anti"
            else: return "planar?"
        except Exception as e:
            return f"Error checking addition config: {e}"

    # --- Main Logic ---
    expected_features = {'ester': 'cis', 'addition': 'anti'}
    
    analysis_results = {}
    errors = []
    for key, smiles in options.items():
        ester_conf = check_ester_config(smiles)
        addition_conf = "N/A"
        
        if "Error" in ester_conf:
            errors.append(f"Option {key}: {ester_conf}")
            analysis_results[key] = {'ester': ester_conf, 'addition': 'N/A'}
            continue
        
        if ester_conf == 'cis':
            addition_conf = check_addition_config(smiles)
            if "Error" in addition_conf:
                errors.append(f"Option {key}: {addition_conf}")
        
        analysis_results[key] = {'ester': ester_conf, 'addition': addition_conf}

    if errors:
        return f"An error occurred during analysis: {'; '.join(errors)}"

    llm_choice_features = analysis_results.get(llm_answer)

    if not llm_choice_features:
        return f"LLM answer '{llm_answer}' is not a valid option."

    # Verify the LLM's choice
    if (llm_choice_features['ester'] == expected_features['ester'] and
        llm_choice_features['addition'] == expected_features['addition']):
        
        # Ensure no other option is also correct
        for key, features in analysis_results.items():
            if key != llm_answer:
                if (features['ester'] == expected_features['ester'] and
                    features['addition'] == expected_features['addition']):
                    return (f"Incorrect. Both option {llm_answer} and option {key} match the expected features. The answer is ambiguous.")
        return "Correct"
    else:
        reasons = []
        if llm_choice_features['ester'] != expected_features['ester']:
            reasons.append(f"it has '{llm_choice_features['ester']}' esters instead of the required '{expected_features['ester']}' configuration")
        if llm_choice_features['addition'] != expected_features['addition']:
            reasons.append(f"it has '{llm_choice_features['addition']}' addition instead of the required '{expected_features['addition']}' configuration")
        
        return (f"Incorrect. The selected answer {llm_answer} is not the major product because " + " and ".join(reasons) + ".")

# Execute the check and print the result.
# The redirect_stdout is used to capture any unwanted prints from libraries.
f = io.StringIO()
with redirect_stdout(f):
    result = check_chemistry_answer()
print(result)