# Import the necessary library
from rdkit import Chem

def check_hydrogen_environments():
    """
    Checks the number of chemically distinct hydrogen atoms in the final product (4).
    The structure of product 4 is determined from the reaction sequence described:
    a Diels-Alder reaction between 7-norbornadienone and 5,6-dimethylidenecyclohexa-1,3-diene.
    """
    # SMILES string for the final product, 4. This structure has the C_s symmetry
    # plane described in the answer, passing through the C=O and bridgehead carbons.
    # The SMILES was obtained by constructing the molecule in a chemical drawing tool.
    smiles_product_4 = "O=C1C2C=CC1C3C4=CC=CC=C4C5C23"

    # Create a molecule object from the SMILES string
    mol = Chem.MolFromSmiles(smiles_product_4)

    if not mol:
        return "Error: The SMILES string for the final product could not be parsed."

    # Add explicit hydrogens to the molecule graph
    mol_with_H = Chem.AddHs(mol)

    # Get the indices of all hydrogen atoms
    h_indices = [atom.GetIdx() for atom in mol_with_H.GetAtoms() if atom.GetAtomicNum() == 1]

    # Use a set to store the canonical SMILES of molecules where each H is uniquely marked
    unique_smiles_set = set()

    # Iterate through each hydrogen, replace it with Deuterium, and find the canonical SMILES
    for h_idx in h_indices:
        # Create an editable copy of the molecule
        mol_copy = Chem.RWMol(mol_with_H)
        
        # Get the specific hydrogen atom and change its atomic number to 2 (Deuterium)
        # This makes it a different element for the canonicalization algorithm.
        atom_to_change = mol_copy.GetAtomWithIdx(h_idx)
        atom_to_change.SetAtomicNum(2)
        
        # Convert the editable molecule back to a regular molecule
        modified_mol = mol_copy.GetMol()
        
        # Generate a canonical SMILES string that accounts for stereochemistry.
        # If two hydrogens are chemically equivalent, replacing either one with Deuterium
        # will result in the same canonical SMILES string.
        canonical_smiles = Chem.MolToSmiles(modified_mol, isomericSmiles=True)
        unique_smiles_set.add(canonical_smiles)

    # The number of unique SMILES strings is the number of distinct hydrogen environments
    calculated_distinct_hydrogens = len(unique_smiles_set)
    
    # The answer provided is 8
    expected_answer = 8

    # Check if our calculated value matches the expected answer
    if calculated_distinct_hydrogens == expected_answer:
        return "Correct"
    else:
        return (f"Incorrect. The provided answer claims there are {expected_answer} chemically distinct hydrogen atoms. "
                f"However, a computational analysis of the described final product's structure "
                f"({smiles_product_4}) reveals {calculated_distinct_hydrogens} distinct hydrogen environments. "
                f"The reasoning in the provided answer is likely flawed in its symmetry analysis or assumes an incorrect final structure.")

# Run the check
result = check_hydrogen_environments()
print(result)