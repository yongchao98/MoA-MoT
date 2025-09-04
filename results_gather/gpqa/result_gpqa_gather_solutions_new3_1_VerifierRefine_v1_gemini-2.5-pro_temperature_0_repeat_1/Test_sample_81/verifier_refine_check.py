# You need to install the RDKit library for this code to run.
# You can install it via pip: pip install rdkit
try:
    from rdkit import Chem
except ImportError:
    print("Error: RDKit is not installed. Please run 'pip install rdkit' to use this checker.")
    exit()

def check_correctness_of_chemistry_answer():
    """
    This function checks the correctness of the provided answer by analyzing the stereochemistry
    of the molecules based on their SMILES strings and comparing it to the reasoning provided in the solution.
    """
    # The SMILES strings for each option as provided in the question.
    smiles_map = {
        'A': "O=C(OC)[C@H]1[C@@H](C=C2)[C@@H]3[C@@H]([C@@H]4[C@H]3[C@@H]5C=C[C@H]4C5)[C@@H]2[C@H]1C(OC)=O",
        'B': "O=C(OC)[C@@H]1[C@@H](C=C2)[C@@H]3[C@@H]([C@@H]4[C@H]3[C@H]5C=C[C@@H]4C5)[C@@H]2[C@H]1C(OC)=O",
        'C': "O=C(OC)[C@@H]1[C@@H](C=C2)[C@@H]3[C@@H]([C@H]4[C@@H]3[C@H]5C=C[C@@H]4C5)[C@@H]2[C@@H]1C(OC)=O",
        'D': "O=C(OC)[C@@H]1[C@H](C=C2)[C@@H]3[C@@H]([C@H]4[C@@H]3[C@@H]5C=C[C@H]4C5)[C@H]2[C@@H]1C(OC)=O"
    }
    
    given_answer = 'A'

    # --- Helper Functions ---

    def get_ester_stereochem(mol):
        """Determines if the two ester groups are cis or trans."""
        # Find the two chiral carbons bonded to the ester groups.
        patt = Chem.MolFromSmarts('[CX4H](C(=O)OC)')
        matches = mol.GetSubstructMatches(patt)
        
        if len(matches) != 2:
            return f"Error: Found {len(matches)} ester-bearing chiral carbons, expected 2."

        atom1_idx, atom2_idx = matches[0][0], matches[1][0]
        tag1 = mol.GetAtomWithIdx(atom1_idx).GetChiralTag()
        tag2 = mol.GetAtomWithIdx(atom2_idx).GetChiralTag()

        # For adjacent atoms in a ring, same chiral tag means cis, different means trans.
        return "cis" if tag1 == tag2 else "trans"

    def get_bridge_orientation(mol):
        """Determines if the bridge orientation is syn or anti based on fusion stereochemistry."""
        Chem.AssignStereochemistry(mol, cleanIt=True, force=True, flagPossibleStereoCenters=True)
        
        # Find the four chiral atoms on the central 4-membered ring.
        fusion_atoms = [a for a in mol.GetAtoms() if a.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED and a.IsInRingSize(4)]

        if len(fusion_atoms) != 4:
            return f"Error: Found {len(fusion_atoms)} fusion atoms on 4-membered ring, expected 4."

        # Find the central 4-membered ring that contains all fusion atoms
        ring4_indices = None
        for r_indices in mol.GetRingInfo().AtomRings():
            if len(r_indices) == 4 and all(a.GetIdx() in r_indices for a in fusion_atoms):
                ring4_indices = list(r_indices)
                break
        if not ring4_indices:
            return "Error: Could not find the central 4-membered ring."

        # Create a map of atom index to R/S label
        label_map = {}
        for atom in fusion_atoms:
            if atom.HasProp('_CIPCode'):
                label_map[atom.GetIdx()] = atom.GetProp('_CIPCode')
            else:
                return "Error: Could not assign R/S labels."
        
        # Order the atoms by their connectivity in the ring
        path = [ring4_indices[0]]
        curr = ring4_indices[0]
        visited = {curr}
        while len(path) < 4:
            found_next = False
            for n in mol.GetAtomWithIdx(curr).GetNeighbors():
                if n.GetIdx() in ring4_indices and n.GetIdx() not in visited:
                    path.append(n.GetIdx())
                    visited.add(n.GetIdx())
                    curr = n.GetIdx()
                    found_next = True
                    break
            if not found_next: return "Error: Could not traverse the 4-membered ring."

        # Count how many adjacent pairs in the ring have the same label
        same_label_neighbors = sum(1 for i in range(4) if label_map[path[i]] == label_map[path[(i + 1) % 4]])
        
        # An alternating pattern (R-S-R-S) has 0 same-label neighbors -> 'anti'.
        # A paired pattern (R-R-S-S) has 2 same-label neighbors -> 'syn'.
        if same_label_neighbors == 0: return "anti"
        elif same_label_neighbors == 2: return "syn"
        else: return "Ambiguous"

    # --- Main Analysis ---
    structural_data = {}
    for key, sm in smiles_map.items():
        mol = Chem.MolFromSmiles(sm)
        if not mol: return f"Error: Invalid SMILES string for option {key}."
        structural_data[key] = {
            'ester': get_ester_stereochem(mol),
            'bridge': get_bridge_orientation(mol)
        }

    # --- Verification Step ---
    # Principles: Major product must have cis-esters and an anti-bridge orientation.
    
    # Check the reasoning in the provided answer against our findings.
    if structural_data['B']['ester'] != 'trans' or structural_data['D']['ester'] != 'trans':
        return "Reasoning check failed: The answer's logic to eliminate B and D because they are 'trans' is not fully supported by the SMILES strings. (Note: B is trans, but D is also trans, so the logic holds)."
    if structural_data['A']['ester'] != 'cis' or structural_data['A']['bridge'] != 'anti':
        return f"Reasoning check failed: The answer identifies A as the cis/anti product, but the code finds it to be {structural_data['A']['ester']}/{structural_data['A']['bridge']}."
    if structural_data['C']['ester'] != 'cis' or structural_data['C']['bridge'] != 'syn':
        return f"Reasoning check failed: The answer identifies C as the cis/syn product, but the code finds it to be {structural_data['C']['ester']}/{structural_data['C']['bridge']}."

    # Determine the correct answer based on chemical principles.
    major_product_key = None
    for key, data in structural_data.items():
        if data['ester'] == 'cis' and data['bridge'] == 'anti':
            major_product_key = key
            
    if major_product_key is None:
        return "Analysis Error: Could not find a structure matching the major product criteria (cis-ester, anti-bridge)."

    if major_product_key == given_answer:
        return "Correct"
    else:
        return f"Incorrect. The analysis shows the major product (cis-ester, anti-bridge) is option {major_product_key}, but the provided answer is {given_answer}."

# Run the check and print the result.
print(check_correctness_of_chemistry_answer())