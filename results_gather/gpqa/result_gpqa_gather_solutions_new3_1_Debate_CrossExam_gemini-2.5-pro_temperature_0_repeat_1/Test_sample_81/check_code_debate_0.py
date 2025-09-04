import sys
from rdkit import Chem
from rdkit.Chem import AllChem
import re

# Suppress RDKit warnings for cleaner output
from rdkit import rdBase
rdBase.DisableLog('rdApp.*')

def get_ester_stereochem(smiles: str) -> str:
    """
    Determines if the two ester groups are cis or trans based on the SMILES string.
    This is a simplified parsing method that works for the given SMILES format.
    The two ester-bearing carbons are denoted by the atom map number 1.
    """
    # Find the stereochemistry of the two carbons attached to the esters.
    # These are marked with the number '1' in the SMILES strings provided.
    matches = re.findall(r'\[C(@{1,2})H\]1', smiles)
    if len(matches) != 2:
        return "Error: Could not identify two ester-bearing carbons."
    
    # If both are @@ or both are @, they are cis. If one is @ and the other is @@, they are trans.
    if matches[0] == matches[1]:
        return "cis"
    else:
        return "trans"

def get_bridge_orientation(smiles: str) -> str:
    """
    Determines if the two main bridges are syn or anti using 3D geometry.
    'anti' means the bridges are on opposite sides of the central ring.
    'syn' means they are on the same side.
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return "Error: Invalid SMILES"

        mol = Chem.AddHs(mol)
        # Embed molecule in 3D space and optimize geometry
        if AllChem.EmbedMolecule(mol, randomSeed=42) == -1:
            return "Error: Could not generate 3D conformer"
        AllChem.MMFFOptimizeMolecule(mol)
        conf = mol.GetConformer()

        # Find the central cyclobutane ring
        patt_c4 = Chem.MolFromSmarts('[#6r4]1[#6r4][#6r4][#6r4]1')
        c4_atoms_match = mol.GetSubstructMatch(patt_c4)
        if not c4_atoms_match:
            return "Error: Could not find central C4 ring"
        c4_atoms = set(c4_atoms_match)

        # Find the methano bridge atom (CH2 from cyclopentadiene)
        patt_methano = Chem.MolFromSmarts('[CH2X4]([#6r4])([#6r4])')
        methano_match = mol.GetSubstructMatch(patt_methano)
        if not methano_match:
            return "Error: Could not find methano bridge"
        methano_C_idx = methano_match[0]
        methano_neighbors = set(mol.GetAtomWithIdx(methano_C_idx).GetNeighbors())
        methano_neighbor_indices = {a.GetIdx() for a in methano_neighbors}

        # The ethano bridge carbons are the other two carbons in the C4 ring
        ethano_bridge_indices = [idx for idx in c4_atoms if idx not in methano_neighbor_indices]
        if len(ethano_bridge_indices) != 2:
            return "Error: Could not identify ethano bridge carbons"
        ethano_C1_idx, ethano_C2_idx = ethano_bridge_indices

        # Calculate the centroid of the C4 ring
        c4_pos = [conf.GetAtomPosition(i) for i in c4_atoms]
        centroid = sum(c4_pos, Chem.rdGeometry.Point3D(0, 0, 0)) / 4.0

        # Get position of the methano bridge atom
        pos_methano = conf.GetAtomPosition(methano_C_idx)
        
        # Get midpoint of the ethano bridge carbons
        pos_ethano_C1 = conf.GetAtomPosition(ethano_C1_idx)
        pos_ethano_C2 = conf.GetAtomPosition(ethano_C2_idx)
        pos_ethano_mid = (pos_ethano_C1 + pos_ethano_C2) / 2.0

        # Calculate vectors from the centroid to the bridges
        v_methano = pos_methano - centroid
        v_ethano = pos_ethano_mid - centroid

        # The dot product determines the relative orientation
        dot_product = v_methano.DotProduct(v_ethano)

        return "anti" if dot_product < 0 else "syn"

    except Exception as e:
        return f"Error: {e}"

def check_answer():
    """
    Main function to check the correctness of the answer.
    """
    options = {
        'A': {
            "name": "dimethyl (1S,4R,4aR,4bR,5S,8R,8aS,8bS,10S,11R)-1,4,4a,4b,5,8,8a,8b-octahydro-1,4-ethano-5,8-methanobiphenylene-10,11-dicarboxylate",
            "smiles": "O=C(OC)[C@@H]1[C@@H](C=C2)[C@@H]3[C@@H]([C@H]4[C@@H]3[C@H]5C=C[C@@H]4C5)[C@@H]2[C@@H]1C(OC)=O"
        },
        'B': {
            "name": "dimethyl (1R,4S,4aS,4bR,5R,8S,8aS,8bR,10R,11R)-1,4,4a,4b,5,8,8a,8b-octahydro-1,4-ethano-5,8-methanobiphenylene-10,11-dicarboxylate",
            "smiles": "O=C(OC)[C@@H]1[C@@H](C=C2)[C@@H]3[C@@H]([C@@H]4[C@H]3[C@H]5C=C[C@@H]4C5)[C@@H]2[C@H]1C(OC)=O"
        },
        'C': {
            "name": "dimethyl (1R,4S,4aR,4bR,5R,8S,8aS,8bS,10S,11R)-1,4,4a,4b,5,8,8a,8b-octahydro-1,4-ethano-5,8-methanobiphenylene-10,11-dicarboxylate",
            "smiles": "O=C(OC)[C@@H]1[C@H](C=C2)[C@@H]3[C@@H]([C@H]4[C@@H]3[C@@H]5C=C[C@H]4C5)[C@H]2[C@@H]1C(OC)=O"
        },
        'D': {
            "name": "dimethyl (1R,4S,4aS,4bR,5S,8R,8aS,8bR,10S,11R)-1,4,4a,4b,5,8,8a,8b-octahydro-1,4-ethano-5,8-methanobiphenylene-10,11-dicarboxylate",
            "smiles": "O=C(OC)[C@H]1[C@@H](C=C2)[C@@H]3[C@@H]([C@@H]4[C@H]3[C@@H]5C=C[C@H]4C5)[C@@H]2[C@H]1C(OC)=O"
        }
    }

    results = []
    correct_candidates = []

    for key, value in options.items():
        ester_config = get_ester_stereochem(value["smiles"])
        bridge_config = get_bridge_orientation(value["smiles"])
        
        is_correct = ester_config == "cis" and bridge_config == "anti"
        results.append(f"Option {key}:")
        results.append(f"  - Ester stereochemistry: {ester_config.upper()}")
        results.append(f"  - Bridge orientation: {bridge_config.upper()}")
        
        if ester_config != "cis":
            results.append("  - Verdict: Incorrect (violates cis-ester constraint from maleic anhydride).")
        elif bridge_config != "anti":
            results.append("  - Verdict: Incorrect (not the major anti-adduct from sterically-favored reaction).")
        else:
            results.append("  - Verdict: Correct (satisfies both stereochemical constraints for the major product).")
            correct_candidates.append(key)
        results.append("-" * 20)

    print("--- Analysis Results ---")
    print("\n".join(results))

    # Final check based on the analysis
    if 'A' in correct_candidates and len(correct_candidates) == 1:
        return "Correct"
    elif 'A' not in correct_candidates:
        return "Incorrect. The provided answer 'A' does not satisfy the stereochemical constraints for the major product according to the analysis."
    else:
        return "Incorrect. The analysis found multiple candidates that satisfy the constraints, indicating a potential ambiguity in the question or options."

# Run the check
result = check_answer()
print(f"\nFinal Conclusion: The provided answer 'A' is {result}")
