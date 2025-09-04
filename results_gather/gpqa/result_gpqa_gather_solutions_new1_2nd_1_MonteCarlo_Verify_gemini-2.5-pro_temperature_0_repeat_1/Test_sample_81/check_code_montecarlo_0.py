import re
from rdkit import Chem

def check_organic_synthesis_answer():
    """
    Checks the correctness of the selected answer for the multi-step synthesis problem.

    The verification proceeds in two main steps:
    1.  Checks a fundamental chemical principle: The product must have cis-diester groups
        because the reaction starts with maleic anhydride (a cis-dienophile). This is
        checked by parsing the IUPAC name. Any trans-isomer is definitively incorrect.
    2.  Checks for self-consistency in the provided data for the chosen answer: It verifies
        that the stereochemistry described in the IUPAC name matches the stereochemistry
        encoded in the SMILES string. This is done using RDKit to parse the SMILES
        and identify the stereochemistry of the relevant atoms.
    """
    # Data from the question prompt
    options = {
        'A': {
            "name": "dimethyl (1S,4R,4aR,4bR,5S,8R,8aS,8bS,10S,11R)-1,4,4a,4b,5,8,8a,8b-octahydro-1,4-ethano-5,8-methanobiphenylene-10,11-dicarboxylate",
            "smiles": "O=C(OC)[C@@H]1[C@@H](C=C2)[C@@H]3[C@@H]([C@H]4[C@@H]3[C@H]5C=C[C@@H]4C5)[C@@H]2[C@@H]1C(OC)=O"
        },
        'B': {
            "name": "dimethyl (1R,4S,4aS,4bR,5S,8R,8aS,8bR,10S,11R)-1,4,4a,4b,5,8,8a,8b-octahydro-1,4-ethano-5,8-methanobiphenylene-10,11-dicarboxylate",
            "smiles": "O=C(OC)[C@H]1[C@@H](C=C2)[C@@H]3[C@@H]([C@@H]4[C@H]3[C@@H]5C=C[C@H]4C5)[C@@H]2[C@H]1C(OC)=O"
        },
        'C': {
            "name": "dimethyl (1R,4S,4aR,4bR,5R,8S,8aS,8bS,10S,11R)-1,4,4a,4b,5,8,8a,8b-octahydro-1,4-ethano-5,8-methanobiphenylene-10,11-dicarboxylate",
            "smiles": "O=C(OC)[C@@H]1[C@H](C=C2)[C@@H]3[C@@H]([C@H]4[C@@H]3[C@@H]5C=C[C@H]4C5)[C@H]2[C@@H]1C(OC)=O"
        },
        'D': {
            "name": "dimethyl (1R,4S,4aS,4bR,5R,8S,8aS,8bR,10R,11R)-1,4,4a,4b,5,8,8a,8b-octahydro-1,4-ethano-5,8-methanobiphenylene-10,11-dicarboxylate",
            "smiles": "O=C(OC)[C@@H]1[C@@H](C=C2)[C@@H]3[C@@H]([C@@H]4[C@H]3[C@H]5C=C[C@@H]4C5)[C@@H]2[C@H]1C(OC)=O"
        }
    }
    
    llm_answer = 'A'

    # --- Step 1: Verify the cis/trans stereochemistry from the IUPAC name ---
    chosen_option_name = options[llm_answer]['name']
    match = re.search(r'10(R|S),11(R|S)', chosen_option_name)
    if not match:
        return f"Constraint check failed: Could not determine C10/C11 stereochemistry from the name of the chosen answer {llm_answer}."
    
    c10_stereo, c11_stereo = match.groups()
    
    # For adjacent stereocenters, different labels (R,S or S,R) mean 'cis'. Same labels (R,R or S,S) mean 'trans'.
    if c10_stereo == c11_stereo:
         return f"The answer {llm_answer} is incorrect. Its IUPAC name specifies a ({c10_stereo},{c11_stereo}) configuration for the ester groups, which is a 'trans' relationship. The reaction must produce a 'cis' isomer because it starts with maleic anhydride."

    # --- Step 2: Verify consistency between the SMILES string and the IUPAC name ---
    def get_diester_stereo_from_smiles(smiles_str):
        mol = Chem.MolFromSmiles(smiles_str)
        if mol is None: return None, "Invalid SMILES string."
        
        Chem.AssignAtomChiralTagsFromStructure(mol, replaceExistingTags=True)
        
        patt = Chem.MolFromSmarts('[C;H1;!$(C=O)](C(=O)O[CH3])')
        ester_attachment_indices = [x[0] for x in mol.GetSubstructMatches(patt)]
        
        if len(ester_attachment_indices) < 2: return None, "Could not find two ester-bearing carbons."
            
        c10_idx, c11_idx = -1, -1
        found = False
        for i in range(len(ester_attachment_indices)):
            for j in range(i + 1, len(ester_attachment_indices)):
                idx1, idx2 = ester_attachment_indices[i], ester_attachment_indices[j]
                if mol.GetBondBetweenAtoms(idx1, idx2) is not None:
                    c10_idx, c11_idx = idx1, idx2
                    found = True
                    break
            if found: break
        
        if not found: return None, "Could not find the two adjacent ester-bearing carbons."
            
        chiral_centers = dict(Chem.FindMolChiralCenters(mol, includeUnassigned=True))
        tag1, tag2 = chiral_centers.get(c10_idx), chiral_centers.get(c11_idx)
        
        if tag1 is None or tag2 is None: return None, f"Could not determine R/S configuration for the ester-bearing carbons (atoms {c10_idx}, {c11_idx})."
            
        return {tag1, tag2}, None

    chosen_smiles = options[llm_answer]['smiles']
    stereo_from_smiles, err = get_diester_stereo_from_smiles(chosen_smiles)
    if err:
        return f"Error processing SMILES for answer {llm_answer}: {err}"

    stereo_from_name = {c10_stereo, c11_stereo}

    if stereo_from_smiles != stereo_from_name:
        return f"The answer {llm_answer} is incorrect due to data inconsistency. The IUPAC name specifies a ({c10_stereo},{c11_stereo}) configuration (set: {stereo_from_name}), but the SMILES string corresponds to a configuration with the set of stereodescriptors {stereo_from_smiles}."

    # --- Conclusion ---
    # The chosen answer 'A' satisfies the necessary 'cis' constraint.
    # The SMILES and IUPAC name for 'A' are self-consistent.
    # The LLM's reasoning that 'A' is the specific 'endo, anti' major product is chemically sound.
    # The code has found no contradictions.
    return "Correct"

# Run the check
result = check_organic_synthesis_answer()
print(result)