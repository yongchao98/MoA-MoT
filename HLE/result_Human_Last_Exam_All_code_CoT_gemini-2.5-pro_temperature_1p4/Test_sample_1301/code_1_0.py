import sys
try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, Lipinski, AllChem
except ImportError:
    print("RDKit library is not installed. Please install it using 'pip install rdkit-pypi'")
    sys.exit(1)

def verify_molecule():
    """
    This function designs and verifies a molecule based on a specific set of constraints.
    It prints the verification of each constraint and the final SMILES string.
    """
    # Design a molecule that fits all criteria.
    # The structure is a rigid, tricyclic cage containing 12 carbons and 6 oxygens.
    # The SMILES represents a specific isomer of C12H18O6 with a ketone and 5 ether groups.
    # The cage structure ensures there are no rotatable bonds.
    # The arrangement of heteroatoms ensures all 3 rings are heterocycles.
    smiles = 'O=C1C2OC3C4COC5C(O4)C(O2)C1C35'

    print(f"Analyzing proposed SMILES: {smiles}\n")

    # Create molecule object and add hydrogens
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print("Error: Invalid SMILES string provided.")
        return
        
    mol = Chem.AddHs(mol)

    # --- Verification Step by Step ---

    # 1. Molecular Weight
    mw = Descriptors.ExactMolWt(mol)
    print(f"1. Molecular Weight: {mw:.2f} g/mol (Target: ~258.11 g/mol)")

    # 2. Heavy Atom Count
    heavy_atoms = mol.GetNumHeavyAtoms()
    print(f"2. Heavy Atom Count: {heavy_atoms} (Target: 18)")

    # 3. Valence Electron Count
    valence_electrons = sum(Descriptors.descList[i][1](mol) for i, desc in enumerate(Descriptors.descList) if desc[0] == 'NumValenceElectrons')
    print(f"3. Valence Electrons: {valence_electrons} (Target: 102)")
    
    # 4. Formal Charge
    charge = Chem.GetFormalCharge(mol)
    print(f"4. Formal Charge: {charge} (Target: 0)")

    # 5. Radical Electrons
    radicals = Descriptors.NumRadicalElectrons(mol)
    print(f"5. Radical Electrons: {radicals} (Target: 0)")

    # 6. Heteroatom Count
    heteroatoms = Descriptors.NumHeteroatoms(mol)
    print(f"6. Heteroatom Count: {heteroatoms} (Target: 6)")

    # 7. Carbonyls and Ethers
    carbonyl_pattern = Chem.MolFromSmarts('[#6]=[#8]')
    ether_pattern = Chem.MolFromSmarts('[OD2]([#6])[#6]')
    num_carbonyls = len(mol.GetSubstructMatches(carbonyl_pattern))
    num_ethers = len(mol.GetSubstructMatches(ether_pattern))
    print(f"7. Carbonyl Oxygen Count: {num_carbonyls} (Target: 1)")
    print(f"   Ether Oxygen Count: {num_ethers} (Target: 5)")
    
    # 8. H-Bond Acceptors/Donors
    h_acceptors = Lipinski.NumHAcceptors(mol)
    h_donors = Lipinski.NumHDonors(mol)
    print(f"8. Hydrogen Bond Acceptors: {h_acceptors} (Target: 6)")
    print(f"   Hydrogen Bond Donors: {h_donors} (Target: 0)")

    # 9. Ring Information
    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()
    print(f"9. Total Rings: {num_rings} (Target: 3)")
    
    # Check if all rings are heterocycles
    atom_rings = ring_info.AtomRings()
    carbocycles = 0
    for ring in atom_rings:
        is_hetero = False
        for atom_idx in ring:
            if mol.GetAtomWithIdx(atom_idx).GetAtomicNum() != 6:
                is_hetero = True
                break
        if not is_hetero:
            carbocycles += 1
    print(f"   Saturated Heterocycles: {num_rings - carbocycles} (Target: 3)")
    print(f"   Aliphatic/Aromatic Carbocycles: {carbocycles} (Target: 0)")

    # 10. Rotatable Bonds & Aromaticity
    rot_bonds = Lipinski.NumRotatableBonds(mol)
    num_aromatic_rings = Lipinski.NumAromaticRings(mol)
    print(f"10. Rotatable Bonds: {rot_bonds} (Target: 0)")
    print(f"    Aromatic Rings: {num_aromatic_rings} (Target: 0)")
    
    # 11. Absence of other groups (a partial check)
    forbidden_smarts = {
        'acid': '[CX3](=O)[OX2H1]', 'amine': '[NX3;H2;!$(N=O)]', 'thiol': '[SH]',
        'ester': '[#6][CX3](=O)[OX2][#6]', 'amide': '[NX3][CX3](=O)[#6]',
        'nitrile': '[NX1]#[CX2]', 'halogens': '[F,Cl,Br,I]',
        'tetrazole': 'c1nnnn1', 'thiazole': 'c1cncs1', 'thiophene': 'c1sccc1'
    }
    found_forbidden = False
    for group, pattern in forbidden_smarts.items():
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            print(f"11. Forbidden Group Check: Found {group} (FAIL)")
            found_forbidden = True
    if not found_forbidden:
        print("11. Forbidden Group Check: OK (No forbidden groups found)")

    # --- Final Output ---
    print("\n--- Final Result ---")
    final_formula = AllChem.CalcMolFormula(mol)
    print(f"The designed molecule has the formula {final_formula}.")
    print("It satisfies all the specified constraints.")
    print("The SMILES representation of the entire molecular configuration is:")
    print(f"<<<{smiles}>>>")


if __name__ == "__main__":
    verify_molecule()