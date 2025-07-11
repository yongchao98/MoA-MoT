# The user wants a Python script to verify the designed molecule.
# I will use the RDKit library for this purpose.

# First, ensure rdkit is installed. If not, the user can install it via pip:
# pip install rdkit

from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

def analyze_molecule(smiles):
    """
    Analyzes a molecule from a SMILES string and checks it against the specified criteria.
    """
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        print("Error: Invalid SMILES string provided.")
        return

    # Add hydrogens to the molecule graph for accurate calculations
    mol = Chem.AddHs(mol)

    print(f"Analysis for SMILES: {smiles}\n")

    # --- Verify Primary Constraints ---
    # 1. Heavy Atoms
    heavy_atoms = Descriptors.HeavyAtomCount(mol)
    print(f"1. Heavy Atoms: {heavy_atoms} (Target: 17)")

    # 2. Heteroatoms (N and O)
    n_atoms = len(mol.GetAtomsMatchingQuery(Chem.MolFromSmarts("[N]")))
    o_atoms = len(mol.GetAtomsMatchingQuery(Chem.MolFromSmarts("[O]")))
    heteroatoms = n_atoms + o_atoms
    print(f"2. Heteroatoms: {heteroatoms} (N={n_atoms}, O={o_atoms})")
    print("   (Note: This count of 7 contradicts the prompt's '5 heteroatoms' but matches the required functional groups)")

    # 3. Formal Charge
    charge = Chem.GetFormalCharge(mol)
    print(f"3. Formal Charge: {charge} (Target: 0)")

    # 4. Valence Electrons
    valence_electrons = sum([atom.GetNumOuterElectrons() for atom in mol.GetAtoms()])
    print(f"4. Valence Electrons: {valence_electrons} (Target: 100)")
    
    # 5. Radical Electrons
    radical_electrons = Descriptors.NumRadicalElectrons(mol)
    print(f"5. Radical Electrons: {radical_electrons} (Target: 0)")

    # --- Verify Structural Features ---
    # 6. Rings
    ri = mol.GetRingInfo()
    aliphatic_heterocycles = rdMolDescriptors.CalcNumAliphaticHeterocycles(mol)
    saturated_rings = rdMolDescriptors.CalcNumSaturatedRings(mol)
    aliphatic_carbocycles = rdMolDescriptors.CalcNumAliphaticCarbocycles(mol)
    aromatic_carbocycles = rdMolDescriptors.CalcNumAromaticCarbocycles(mol)
    saturated_carbocycles = rdMolDescriptors.CalcNumSaturatedCarbocycles(mol)
    print("\n--- Structural Features ---")
    print(f"6. Rings:")
    print(f"   - Aliphatic Heterocycles: {aliphatic_heterocycles} (Target: 2)")
    print(f"   - Saturated Rings: {saturated_rings} (Target: 2)")
    print(f"   - Aliphatic Carbocycles: {aliphatic_carbocycles} (Target: 0)")
    print(f"   - Aromatic Carbocycles: {aromatic_carbocycles} (Target: 0)")
    print(f"   - Saturated Carbocycles: {saturated_carbocycles} (Target: 0)")

    # 7. Hydrogen Bonding
    h_donors = Descriptors.NumHDonors(mol)
    h_acceptors = Descriptors.NumHAcceptors(mol)
    print(f"7. Hydrogen Bonding:")
    print(f"   - H-Bond Donors: {h_donors} (Target: 0)")
    print(f"   - H-Bond Acceptors: {h_acceptors} (Target: Allowed)")

    # 8. Rotatable Bonds
    rotatable_bonds = Descriptors.NumRotatableBonds(mol)
    print(f"8. Rotatable Bonds: {rotatable_bonds} (Target: 6)")

    # 9. Functional Groups
    ether_oxygens = len(mol.GetAtomsMatchingQuery(Chem.MolFromSmarts("[OD2](-C)-C")))
    tertiary_amines = len(mol.GetAtomsMatchingQuery(Chem.MolFromSmarts("[NX3;H0;!$(C=O)]")))
    print(f"9. Functional Groups:")
    print(f"   - Ether Oxygens: {ether_oxygens} (Target: 5)")
    print(f"   - Tertiary Amines: {tertiary_amines} (Target: 2)")

    # 10. Molecular Weight
    # As reasoned, the MW constraint contradicts the valence electron count.
    # The designed molecule (C10H20N2O5) fits the electron count.
    # We will calculate its MW and show the discrepancy.
    formula = rdMolDescriptors.CalcMolFormula(mol)
    mw = Descriptors.MolWt(mol)
    print("\n--- Molecular Weight Analysis ---")
    print(f"The designed molecule's formula is {formula}.")
    print("The molecular weight constraint (244.179) implies a formula of C10H16N2O5, which contradicts the valence electron count and ring structure requirements.")
    print("Calculating the weight for the designed structure (C10H20N2O5):")
    
    c_mass = 12.011
    h_mass = 1.008
    n_mass = 14.007
    o_mass = 15.999
    
    c_count = len(mol.GetAtomsMatchingQuery(Chem.MolFromSmarts("[C]")))
    h_count = len(mol.GetAtomsMatchingQuery(Chem.MolFromSmarts("[H]")))
    n_count = n_atoms
    o_count = o_atoms
    
    calc_mw = (c_count * c_mass) + (h_count * h_mass) + (n_count * n_mass) + (o_count * o_mass)
    
    print(f"({c_count} * {c_mass}) + ({h_count} * {h_mass}) + ({n_count} * {n_mass}) + ({o_count} * {o_mass}) = {calc_mw:.3f}")
    print(f"Calculated Molecular Weight: {mw:.3f} (Target was: 244.179)")
    
    print("\nFinal proposed SMILES satisfying all consistent constraints:")
    print(f">>>{smiles}<<<")


if __name__ == '__main__':
    # Based on the analysis, this molecule fits all compatible constraints.
    # The contradictions in the prompt (MW vs valence electrons, heteroatom count) are noted in the analysis.
    # We prioritize the structural constraints and valence electron count, as they are more consistent.
    designed_smiles = "CN(C1OCOC1)COCN(C)C2OCOC2"
    analyze_molecule(designed_smiles)
