# To run this script, you need to have the RDKit library installed.
# You can install it via pip: pip install rdkit-pypi

from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, rdMolDescriptors

def verify_molecule():
    """
    This function designs and verifies a molecule based on a specific set of chemical rules.
    It prints the verification for each rule and finally outputs the molecule's SMILES string.
    """
    # Proposed SMILES string based on the analysis
    smiles = "Oc1ccc(N=C(C)c2cn(CCC)[nH]c2)cc1"
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)

    print("--- Verifying Molecule Properties ---")
    print(f"Proposed SMILES: {smiles}\n")

    # 1. Heavy Atom Count
    heavy_atoms = mol.GetNumHeavyAtoms()
    print(f"1. Heavy Atom Count: {heavy_atoms} (Target: 18)")
    assert heavy_atoms == 18

    # 2. Molecular Weight
    mw = Descriptors.ExactMolWt(mol)
    print(f"2. Molecular Weight: {mw:.5f} (Target: 243.137)")
    assert abs(mw - 243.137) < 0.001

    # 3. Formal Charge
    charge = Chem.GetFormalCharge(mol)
    print(f"3. Formal Charge: {charge} (Target: 0)")
    assert charge == 0

    # 4. Valence Electron Count
    # C14H17N3O
    c_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    h_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 1)
    n_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    o_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    valence_electrons = (c_atoms * 4) + (h_atoms * 1) + (n_atoms * 5) + (o_atoms * 6)
    print(f"4. Valence Electron Count from Formula ({c_atoms}*4 + {h_atoms}*1 + {n_atoms}*5 + {o_atoms}*6): {valence_electrons} (Target: 94)")
    assert valence_electrons == 94
    
    # 5. Rings (2 aromatic: 1 benzene, 1 imidazole)
    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()
    is_benzene = mol.HasSubstructMatch(Chem.MolFromSmarts("c1ccccc1"))
    is_imidazole = mol.HasSubstructMatch(Chem.MolFromSmarts("c1cn[nH]c1"))
    print(f"5. Aromatic Rings: Found {num_rings} rings. Benzene: {is_benzene}, Imidazole: {is_imidazole} (Target: True for both)")
    assert num_rings == 2 and is_benzene and is_imidazole

    # 6. Functional Groups and Heteroatoms
    num_h_donors = Lipinski.NumHDonors(mol)
    # Lipinski acceptors (N+O count) is a common definition
    num_h_acceptors = Lipinski.NumHAcceptors(mol)
    
    phenolic_oh = mol.HasSubstructMatch(Chem.MolFromSmarts("[OH]c1ccccc1"))
    imine = mol.HasSubstructMatch(Chem.MolFromSmarts("[#6X3]=[#7X2]"))
    # The assumption is one tertiary amine, not three
    tertiary_amine = len(mol.GetSubstructMatches(Chem.MolFromSmarts("[#7X3;H0;!$(*=[O,S,N])]")))
    
    print(f"6. Functional Groups & Properties:")
    print(f"   - Phenolic OH: {phenolic_oh} (Target: True)")
    print(f"   - Imine group: {imine} (Target: True)")
    print(f"   - Tertiary amines: {tertiary_amine} (Target: 1, based on interpretation)")
    print(f"   - Hydrogen Bond Donors: {num_h_donors} (Target: 1)")
    print(f"   - Hydrogen Bond Acceptors: {num_h_acceptors} (Target: 4)")
    assert phenolic_oh and imine and tertiary_amine == 1
    assert num_h_donors == 1 and num_h_acceptors == 4

    # 7. Rotatable Bonds
    rotatable_bonds = Lipinski.NumRotatableBonds(mol)
    print(f"7. Rotatable Bonds: {rotatable_bonds} (Target: 5)")
    assert rotatable_bonds == 5
    
    print("\n--- All Checks Passed ---")
    print("\nFinal Proposed Molecule SMILES:")
    print(smiles)


if __name__ == '__main__':
    verify_molecule()
<<<Oc1ccc(N=C(C)c2cn(CCC)[nH]c2)cc1>>>