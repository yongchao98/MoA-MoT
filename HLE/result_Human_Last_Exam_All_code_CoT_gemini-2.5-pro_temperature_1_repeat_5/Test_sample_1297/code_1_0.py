try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, rdMolDescriptors
except ImportError:
    print("RDKit library is not installed. Please install it using 'pip install rdkit'")
    exit()

def analyze_molecule(smiles):
    """
    Analyzes a molecule from a SMILES string and prints its properties based on the user's criteria.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print(f"Error: Invalid SMILES string provided: {smiles}")
        return

    # Add hydrogens to get correct properties for valence electrons, MW, etc.
    mol = Chem.AddHs(mol)

    # --- Calculate Properties ---
    heavy_atoms = mol.GetNumHeavyAtoms()
    heteroatoms = Descriptors.NumHeteroatoms(mol)
    formal_charge = Chem.GetFormalCharge(mol)
    valence_electrons = Descriptors.NumValenceElectrons(mol)
    radical_electrons = Descriptors.NumRadicalElectrons(mol)
    
    # Ring counts
    ring_info = mol.GetRingInfo()
    aliphatic_heterocycles = Descriptors.NumAliphaticHeterocycles(mol)
    saturated_rings = Descriptors.NumSaturatedRings(mol)
    aliphatic_carbocycles = Descriptors.NumAliphaticCarbocycles(mol)
    aromatic_carbocycles = Descriptors.NumAromaticCarbocycles(mol)
    saturated_carbocycles = Descriptors.NumSaturatedCarbocycles(mol)

    # Hydrogen bonding
    h_bond_donors = Descriptors.NumHDonors(mol)
    h_bond_acceptors = Descriptors.NumHAcceptors(mol)
    
    # Other properties
    rotatable_bonds = Descriptors.NumRotatableBonds(mol)
    molecular_weight = Descriptors.MolWt(mol)
    molecular_formula = rdMolDescriptors.CalcMolFormula(mol)

    # Functional group counts using SMARTS patterns
    ether_pattern = Chem.MolFromSmarts('[OD2](-!@[#6])-!@[#6]')
    tertiary_amine_pattern = Chem.MolFromSmarts('[NX3](-!@[#6])(-!@[#6])(-!@[#6])')
    
    ether_oxygens = len(mol.GetSubstructMatches(ether_pattern))
    tertiary_amines = len(mol.GetSubstructMatches(tertiary_amine_pattern))

    # --- Print Report ---
    print("--- Analysis of Designed Molecule ---")
    print(f"Proposed SMILES: {smiles}")
    print(f"Molecular Formula: {molecular_formula}")
    print("\n--- Verifying Constraints ---")
    print(f"Total Heavy Atoms: {heavy_atoms} (Target: 17)")
    print(f"Total Heteroatoms: {heteroatoms} (*Note: Target of 5 contradicts functional group needs of 7)")
    print(f"Formal Charge: {formal_charge} (Target: 0)")
    print(f"Valence Electrons: {valence_electrons} (Target: 100)")
    print(f"Radical Electrons: {radical_electrons} (Target: 0)")
    print(f"Molecular Weight: {molecular_weight:.3f} (*Note: Target of 244.179 contradicts electron count)")
    
    print("\n--- Ring Structure ---")
    print(f"Aliphatic Heterocycles: {aliphatic_heterocycles} (Target: 2)")
    print(f"Saturated Rings: {saturated_rings} (Target: 2)")
    print(f"Aliphatic Carbocycles: {aliphatic_carbocycles} (Target: 0)")
    print(f"Saturated Carbocycles: {saturated_carbocycles} (Target: 0)")
    
    print("\n--- Hydrogen Bonding ---")
    print(f"Hydrogen Bond Donors: {h_bond_donors} (Target: 0)")
    print(f"Hydrogen Bond Acceptors: {h_bond_acceptors} (Target: available)")

    print("\n--- Functional Groups & Bonds ---")
    print(f"Rotatable Bonds: {rotatable_bonds} (Target: 6)")
    print(f"Ether Oxygens: {ether_oxygens} (Target: 5)")
    print(f"Tertiary Amines: {tertiary_amines} (Target: 2)")
    print("\n--- Final Proposed SMILES ---")

# The SMILES string for the designed molecule:
# (1,4-dioxan-2-yl)-N(CH3)-CH2-N(CH3)-(1,3,5-trioxan-2-yl)
# A more canonical representation might be generated, but this is chemically correct.
# A canonical SMILES for this is: CN(CC1COCCO1)CN(C)C2OCOCO2
proposed_smiles = "CN(C1OCCOC1)CN(C2OCOCO2)C"

# Run the analysis
analyze_molecule(proposed_smiles)

# The final answer as requested
print(f"\n<<<{proposed_smiles}>>>")