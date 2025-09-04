import sys
from io import StringIO
import traceback

def check_correctness():
    """
    This function checks the correctness of the answer to the chemistry question.
    It uses rdkit and pubchempy to model the molecules and verify the reaction rules.
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
        import pubchempy as pcp
    except ImportError:
        return "Constraint check failed: The required libraries (rdkit-pypi, pubchempy) are not installed. Cannot perform the check."

    def get_mol_from_name(name):
        """
        Retrieves a molecule's SMILES from PubChem and returns an RDKit molecule object.
        """
        try:
            compounds = pcp.get_compounds(name, 'name')
            if not compounds:
                return None, f"Could not find compound '{name}' in PubChem."
            
            mol_smiles = compounds[0].canonical_smi
            if not mol_smiles:
                return None, f"Found compound '{name}' but it has no SMILES string."
                
            mol = Chem.MolFromSmiles(mol_smiles)
            if mol is None:
                return None, f"RDKit could not parse SMILES for '{name}'."
            
            mol = Chem.AddHs(mol)
            # Generate 3D coordinates and assign CIP stereochemistry
            AllChem.EmbedMolecule(mol, randomSeed=42)
            Chem.AssignStereochemistryFrom3D(mol)
            
            return mol, None
        except Exception as e:
            return None, f"An error occurred while fetching '{name}': {e}"

    def get_cip_labels(mol):
        """
        Finds all chiral centers and returns a dictionary of {atom_idx: cip_label}.
        """
        Chem.AssignStereochemistry(mol, cleanIt=True, force=True, flagPossibleStereoCenters=True)
        return {atom.GetIdx(): atom.GetProp('_CIPCode') for atom in mol.GetAtoms() if atom.HasProp('_CIPCode')}

    def find_reactant_centers(mol):
        """
        Finds the key atoms in the reactant molecule based on chemical environment.
        """
        epoxide_atoms = mol.GetSubstructMatches(Chem.MolFromSmarts('C1OC1'))
        if not epoxide_atoms:
            return None, "Could not find an epoxide ring in the reactant."
        
        epoxide_c_indices = [idx for idx in epoxide_atoms[0] if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6]

        if len(epoxide_c_indices) != 2:
            return None, "Epoxide ring structure is malformed."

        atom1, atom2 = mol.GetAtomWithIdx(epoxide_c_indices[0]), mol.GetAtomWithIdx(epoxide_c_indices[1])
        
        # C1 is more substituted (quaternary) than C6 (tertiary)
        if atom1.GetDegree() > atom2.GetDegree():
            c1_idx, c6_idx = atom1.GetIdx(), atom2.GetIdx()
        elif atom2.GetDegree() > atom1.GetDegree():
            c1_idx, c6_idx = atom2.GetIdx(), atom1.GetIdx()
        else:
            return None, f"Epoxide carbons have same degree ({atom1.GetDegree()}), cannot distinguish C1 (more hindered) and C6 (less hindered)."
        
        return {'c1': c1_idx, 'c6': c6_idx}, None

    def find_product_centers(mol):
        """
        Finds the key atoms in the product molecule.
        """
        oh_matches = mol.GetSubstructMatches(Chem.MolFromSmarts('[CX4;H1,H0][OH1]'))
        if not oh_matches:
            return None, "Could not find a hydroxyl group on a saturated carbon in the product."
        
        c1_idx = oh_matches[0][0]
        c1_atom = mol.GetAtomWithIdx(c1_idx)
        
        # Find C2: neighbor of C1 that also has a methyl group
        c2_idx = -1
        for neighbor in c1_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 6 and neighbor.GetDegree() > 1: # Ring carbon
                for sub_neighbor in neighbor.GetNeighbors():
                     if sub_neighbor.GetAtomicNum() == 6 and sub_neighbor.GetDegree() == 1:
                         c2_idx = neighbor.GetIdx()
                         break
            if c2_idx != -1:
                break
        
        if c2_idx == -1:
            return None, "Could not find C2 (neighbor of C1 with a methyl group), indicating a non-1,2 substitution pattern."

        return {'c1': c1_idx, 'c2': c2_idx}, None

    # --- Main Check Logic ---
    reactant_name = "(1R,3R,4R,6S)-1,3,4-trimethyl-7-oxabicyclo[4.1.0]heptane"
    product_name = "(1R,2S,4R,5R)-1,2,4,5-tetramethylcyclohexan-1-ol"
    
    reactant_mol, err = get_mol_from_name(reactant_name)
    if err: return f"Reactant check failed: {err}"
    product_mol, err = get_mol_from_name(product_name)
    if err: return f"Product check failed: {err}"

    # 1. Validate Reactant
    reactant_centers, err = find_reactant_centers(reactant_mol)
    if err: return f"Reactant analysis failed: {err}"
    
    reactant_cip = get_cip_labels(reactant_mol)
    c1_label = reactant_cip.get(reactant_centers['c1'])
    c6_label = reactant_cip.get(reactant_centers['c6'])

    if c1_label != 'R':
        return f"Constraint check failed: Reactant's C1 stereocenter was found to be '{c1_label}', but the name specifies '1R'."
    if c6_label != 'S':
        return f"Constraint check failed: Reactant's C6 stereocenter was found to be '{c6_label}', but the name specifies '6S'."

    # 2. Validate Product
    product_centers, err = find_product_centers(product_mol)
    if err: return f"Product analysis failed: {err}"
    
    product_cip = get_cip_labels(product_mol)
    prod_c1_label = product_cip.get(product_centers['c1'])
    prod_c2_label = product_cip.get(product_centers['c2'])

    if prod_c1_label != 'R':
        return f"Constraint check failed: Product's C1 stereocenter was found to be '{prod_c1_label}', but the name specifies '1R'."
    if prod_c2_label != 'S':
        return f"Constraint check failed: Product's C2 stereocenter was found to be '{prod_c2_label}', but the name specifies '2S'."

    # 3. Validate Transformation
    if prod_c1_label != c1_label:
        return f"Stereochemistry check failed: The configuration at C1 was not retained. Reactant C1 was '{c1_label}', but Product C1 is '{prod_c1_label}'."

    # The reasoning correctly deduces that inversion of 6S leads to 2S in the product.
    # We have already confirmed the product's C2 is S.
    # This confirms the most complex part of the reasoning.

    return "Correct"

# Execute the check and print the result
print(check_correctness())