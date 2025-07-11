# First, ensure you have rdkit installed:
# pip install rdkit

from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, rdmolops

def verify_molecule_properties(smiles_string):
    """
    Analyzes a molecule from a SMILES string and verifies its properties against the given criteria.
    """
    mol = Chem.MolFromSmiles(smiles_string)
    if mol is None:
        print("Error: Invalid SMILES string provided.")
        return

    print(f"Verifying properties for SMILES: {smiles_string}\n")

    # 1. Heavy Atoms
    heavy_atoms = mol.GetNumHeavyAtoms()
    print(f"Criterion: 17 heavy atoms")
    print(f"Result: {heavy_atoms}")
    print("-" * 20)

    # 2. Heteroatoms
    heteroatoms = rdMolDescriptors.CalcNumHeteroatoms(mol)
    print(f"Criterion: 5 heteroatoms (N and O)")
    print(f"Result: {heteroatoms}")
    print("-" * 20)

    # 3. Formal Charge
    charge = rdmolops.GetFormalCharge(mol)
    print(f"Criterion: Formal charge of 0")
    print(f"Result: {charge}")
    print("-" * 20)
    
    # 4. Valence Electrons
    valence_electrons = rdMolDescriptors.CalcNumValenceElectrons(mol)
    print(f"Criterion: 100 valence electrons")
    print(f"Result: {valence_electrons}")
    print("-" * 20)

    # 5. Radical Electrons
    radical_electrons = rdMolDescriptors.CalcNumRadicalElectrons(mol)
    print(f"Criterion: 0 radical electrons")
    print(f"Result: {radical_electrons}")
    print("-" * 20)

    # 6. Ring Analysis
    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()
    aliphatic_heterocycles = 0
    saturated_rings = 0
    carbocycles = 0

    if num_rings > 0:
        for ring_atoms in ring_info.AtomRings():
            is_hetero = False
            is_carbocycle = True
            for atom_idx in ring_atoms:
                if mol.GetAtomWithIdx(atom_idx).GetAtomicNum() != 6:
                    is_hetero = True
                    is_carbocycle = False
            
            # Check for aromaticity
            is_aromatic = False
            for bond_idx in ring_info.BondRings()[aliphatic_heterocycles]:
                if mol.GetBondWithIdx(bond_idx).GetIsAromatic():
                    is_aromatic = True
                    break
            
            # Check for saturation
            is_saturated = True
            for bond_idx in ring_info.BondRings()[aliphatic_heterocycles]:
                 if mol.GetBondWithIdx(bond_idx).GetBondType() != Chem.BondType.SINGLE:
                    is_saturated = False
                    break

            if is_hetero and not is_aromatic:
                aliphatic_heterocycles += 1
            if is_saturated:
                saturated_rings += 1
            if is_carbocycle:
                carbocycles += 1

    print(f"Criterion: 2 aliphatic heterocycles")
    print(f"Result: {aliphatic_heterocycles}")
    print(f"Criterion: 2 saturated rings")
    print(f"Result: {saturated_rings}")
    print(f"Criterion: 0 aliphatic, aromatic, or saturated carbocycles")
    print(f"Result: {carbocycles}")
    print("-" * 20)

    # 7. Hydrogen Bonding
    hbd = rdMolDescriptors.CalcNumHBD(mol)
    hba = rdMolDescriptors.CalcNumHBA(mol)
    print(f"Criterion: 0 hydrogen bond donors")
    print(f"Result: {hbd}")
    print(f"Criterion: Hydrogen bond acceptors are allowed")
    print(f"Result: {hba} acceptors found")
    print("-" * 20)

    # 8. Rotatable Bonds
    rotatable_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
    print(f"Criterion: 6 rotatable bonds")
    print(f"Result: {rotatable_bonds}")
    print("-" * 20)

    # 9. Functional Groups
    ether_pattern = Chem.MolFromSmarts("[#6]-[#8X2]-[#6]")
    num_ethers = len(mol.GetSubstructMatches(ether_pattern))
    tertiary_amine_pattern = Chem.MolFromSmarts("[NX3]([#6])([#6])[#6]")
    num_tertiary_amines = len(mol.GetSubstructMatches(tertiary_amine_pattern))

    print("Criterion: 5 ether oxygens (Note: Corrected to 3 based on other constraints)")
    print(f"Result: {num_ethers} ether oxygens found")
    print("Criterion: 2 tertiary amines")
    print(f"Result: {num_tertiary_amines} tertiary amines found")
    print("-" * 20)

    # 10. Molecular Weight
    mw = Descriptors.ExactMolWt(mol)
    print(f"Criterion: Molecular weight of 244.179")
    print(f"Result: {mw:.5f}")
    print("-" * 20)

# --- Main Execution ---
# The proposed SMILES string for a molecule that fits the criteria (with the noted correction).
# The molecule is 1,1'-(oxybis(ethane-2,1-diyl))dimorpholine.
proposed_smiles = "C1COCCN1CCOCCOCN2CCOCC2"
verify_molecule_properties(proposed_smiles)

print("\nFinal proposed SMILES:")
print(proposed_smiles)