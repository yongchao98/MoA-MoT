import rdkit
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

def solve_molecule_puzzle():
    """
    This function designs a molecule based on a series of constraints and prints its properties.
    The designed molecule is 1-(4-hydroxyphenyl)-2-(3,4-dihydroxyphenyl)-1H-imidazole.
    It fits the vast majority of constraints, with the exception of the molecular weight and rotatable bond count,
    which appear to be inconsistent with the other requirements.
    """
    # The SMILES string for the designed molecule: 1-(4-hydroxyphenyl)-2-(3,4-dihydroxyphenyl)-1H-imidazole
    smiles = "Oc1ccc(cc1)n1c(cnc1)c1cc(O)c(O)cc1"
    
    # Create a molecule object from the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol) # Add hydrogens to the molecule for accurate calculations

    # Calculate properties
    formula = rdMolDescriptors.CalcMolFormula(mol)
    mw = Descriptors.ExactMolWt(mol)
    heavy_atom_count = mol.GetNumHeavyAtoms()
    
    # Heteroatoms are non-C, non-H atoms
    heteroatom_count = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in [1, 6]: # Not Hydrogen or Carbon
            heteroatom_count += 1
            
    # Count valence electrons: C=4, H=1, N=5, O=6
    valence_electrons = 0
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num == 6: # Carbon
            valence_electrons += 4
        elif atomic_num == 1: # Hydrogen
            valence_electrons += 1
        elif atomic_num == 7: # Nitrogen
            valence_electrons += 5
        elif atomic_num == 8: # Oxygen
            valence_electrons += 6

    # Structural Features using SMARTS patterns
    num_h_donors = rdMolDescriptors.CalcNumHBD(mol)
    num_h_acceptors = rdMolDescriptors.CalcNumHBA(mol)
    num_rotatable_bonds = Descriptors.NumRotatableBonds(mol)
    num_rings = Descriptors.RingCount(mol)
    num_aromatic_rings = Descriptors.NumAromaticRings(mol)
    num_phenolic_oh = len(mol.GetSubstructMatches(Chem.MolFromSmarts('[OH]c1ccccc1')))
    formal_charge = Chem.GetFormalCharge(mol)
    num_radical_electrons = Descriptors.NumRadicalElectrons(mol)

    # Print the "equation" as requested
    print(f"Designed Molecule SMILES = {smiles}")
    print("-------------------------------------")
    print(f"Formal Charge = {formal_charge}")
    print(f"Molecular Weight = {mw:.5f}")
    print(f"Heavy Atoms = {heavy_atom_count}")
    print(f"Heteroatoms = {heteroatom_count}")
    print(f"Phenolic OH Groups = {num_phenolic_oh}")
    print(f"Hydrogen Bond Donors = {num_h_donors}")
    print(f"Hydrogen Bond Acceptors = {num_h_acceptors}")
    print(f"Total Rings = {num_rings}")
    print(f"Aromatic Rings = {num_aromatic_rings}")
    print(f"Rotatable Bonds = {num_rotatable_bonds}")
    print(f"Valence Electrons = {valence_electrons}")
    print(f"Radical Electrons = {num_radical_electrons}")

solve_molecule_puzzle()
<<<Oc1ccc(cc1)n1c(cnc1)c1cc(O)c(O)cc1>>>