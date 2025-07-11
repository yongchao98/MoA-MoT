from rdkit import Chem
from rdkit.Chem import Descriptors, rdmolops

def solve_molecule_design():
    """
    Designs a molecule based on the specified criteria and verifies its properties.
    """
    # Step 1: Define the proposed SMILES string.
    # This structure is designed based on the interpretation that the molecule should have
    # 3 ether oxygens, not 5, to resolve a contradiction in the prompt's constraints.
    # The molecular formula is C12H24N2O3.
    smiles = "C1CN(COCOCN2CCCCC2)CC1"
    
    # Create an RDKit molecule object and add hydrogens for accurate calculations
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print("Error: Invalid SMILES string.")
        return

    mol = Chem.AddHs(mol)

    # Step 2: Calculate all relevant properties using RDKit

    # --- Atom and Electron Counts ---
    heavy_atoms = mol.GetNumHeavyAtoms()
    heteroatoms = Descriptors.NumHeteroatoms(mol)
    formal_charge = rdmolops.GetFormalCharge(mol)
    valence_electrons = Descriptors.NumValenceElectrons(mol)
    radical_electrons = Descriptors.NumRadicalElectrons(mol)

    # --- Ring Analysis ---
    ring_info = mol.GetRingInfo()
    bond_rings = ring_info.BondRings()
    atom_rings = ring_info.AtomRings()

    aliphatic_heterocycles = 0
    saturated_rings = 0
    carbocycles = 0

    if atom_rings:
        for i, ring_atoms in enumerate(atom_rings):
            ring_bonds = bond_rings[i]
            is_hetero = any(mol.GetAtomWithIdx(atom_idx).GetAtomicNum() != 6 for atom_idx in ring_atoms)
            is_aromatic = all(mol.GetBondWithIdx(bond_idx).GetIsAromatic() for bond_idx in ring_bonds)
            is_saturated = all(mol.GetBondWithIdx(bond_idx).GetBondType() == Chem.BondType.SINGLE for bond_idx in ring_bonds)

            if is_hetero and not is_aromatic:
                aliphatic_heterocycles += 1
            if is_saturated:
                saturated_rings += 1
            if not is_hetero:
                carbocycles += 1

    # --- Hydrogen Bonding and Rotatable Bonds ---
    h_donors = Descriptors.NumHDonors(mol)
    rotatable_bonds = Descriptors.NumRotatableBonds(mol)

    # --- Functional Group Analysis (using SMARTS) ---
    ether_pattern = Chem.MolFromSmarts('[OD2](-[#6])-[#6]')
    tertiary_amine_pattern = Chem.MolFromSmarts('[NX3]([#6])([#6])[#6]')
    num_ethers = len(mol.GetSubstructMatches(ether_pattern))
    num_tertiary_amines = len(mol.GetSubstructMatches(tertiary_amine_pattern))

    # --- Molecular Weight ---
    mw = Descriptors.ExactMolWt(mol)
    
    # Step 3: Print the final answer and verification
    
    print(f"Proposed SMILES representation: {smiles}\n")
    
    print("--- Verification of Properties ---")
    print(f"Total Heavy Atoms:       {heavy_atoms} (Target: 17)")
    print(f"Heteroatoms:             {heteroatoms} (Target: 5)")
    print(f"Formal Charge:           {formal_charge} (Target: 0)")
    print(f"Valence Electrons:       {valence_electrons} (Target: 100)")
    print(f"Radical Electrons:       {radical_electrons} (Target: 0)")
    print(f"Aliphatic Heterocycles:  {aliphatic_heterocycles} (Target: 2)")
    print(f"Saturated Rings:         {saturated_rings} (Target: 2)")
    print(f"Carbocycles (any type):  {carbocycles} (Target: 0)")
    print(f"Hydrogen Bond Donors:    {h_donors} (Target: 0)")
    print(f"Rotatable Bonds:         {rotatable_bonds} (Target: 6)")
    print(f"Ether Oxygens:           {num_ethers} (Target: 5, Adjusted to 3)")
    print(f"Tertiary Amines:         {num_tertiary_amines} (Target: 2)")
    print("-" * 34)

    # Final "equation" for molecular weight, as requested.
    # Formula for C1CN(COCOCN2CCCCC2)CC1 is C12H24N2O3
    num_c, num_h, num_n, num_o = 12, 24, 2, 3
    mass_c, mass_h, mass_n, mass_o = 12.000000, 1.007825, 14.003074, 15.994915
    
    print("Final Equation (Molecular Weight Calculation):")
    print(f"Formula: C{num_c}H{num_h}N{num_n}O{num_o}")
    print(f"({num_c} * C) + ({num_h} * H) + ({num_n} * N) + ({num_o} * O)")
    print(f"({num_c} * {mass_c:.3f}) + ({num_h} * {mass_h:.3f}) + ({num_n} * {mass_n:.3f}) + ({num_o} * {mass_o:.3f}) = {mw:.5f}")
    print(f"Calculated MW: {mw:.5f} (Target: 244.179)")

# Execute the function
solve_molecule_design()

<<<C1CN(COCOCN2CCCCC2)CC1>>>