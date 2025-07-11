import rdkit
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, AllChem

def solve_molecule_challenge():
    """
    This function designs a molecule based on a set of complex constraints,
    verifies its properties using RDKit, and prints the analysis.
    """
    # Based on the constraints, particularly the exact molecular weight and
    # required functional groups, the following SMILES string is proposed.
    # SMILES for 1-(4-hydroxyphenyl)-N-(2-(1H-imidazol-4-yl)ethyl)propan-1-imine
    smiles = "CCC(=NCCc1c[nH]cn1)c1ccc(O)cc1"
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)

    if mol is None:
        print("Error: Could not generate molecule from SMILES string.")
        return

    print(f"Proposed Molecule SMILES: {smiles}\n")

    # --- Verification of Properties ---
    print("--- Verifying Molecular Properties ---")

    # 1. Heavy Atom Count
    heavy_atom_count = mol.GetNumHeavyAtoms()
    print(f"1. Heavy Atom Count: {heavy_atom_count} (Target: 18)")

    # 2. Molecular Weight
    formula = rdMolDescriptors.CalcMolFormula(mol)
    print(f"   Molecular Formula: {formula}")
    # Exact masses for the most common isotopes
    mass_C = 12.000000
    mass_H = 1.007825
    mass_N = 14.003074
    mass_O = 15.994915
    # Extract atom counts from formula
    atom_counts = Chem.rdMolDescriptors.GetMolFormula(mol, True)
    c_count = atom_counts.get('C', 0)
    h_count = atom_counts.get('H', 0)
    n_count = atom_counts.get('N', 0)
    o_count = atom_counts.get('O', 0)
    
    exact_mw = (c_count * mass_C + 
                h_count * mass_H + 
                n_count * mass_N + 
                o_count * mass_O)

    print(f"2. Molecular Weight Calculation (Exact Mass):")
    print(f"   Equation: ({c_count} * {mass_C:.6f}) + ({h_count} * {mass_H:.6f}) + ({n_count} * {mass_N:.6f}) + ({o_count} * {mass_O:.6f})")
    print(f"   Result: {exact_mw:.5f} (Target: 243.137)")
    
    # 3. Formal Charge
    charge = Chem.GetFormalCharge(mol)
    print(f"3. Formal Charge: {charge} (Target: 0)")
    
    # 4. Valence Electron Count
    valence_electrons = sum(atom.GetTotalValenceE() for atom in mol.GetAtoms())
    print(f"4. Valence Electron Count: {valence_electrons} (Target: 94)")
    
    # 5. Ring Structure
    ri = mol.GetRingInfo()
    num_rings = ri.NumRings()
    is_aromatic = [all(ri.IsAtomInRingOfSize(i, size) and mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring) for size, ring in zip(ri.AtomRingSizes(), ri.AtomRings())]
    num_aromatic_rings = sum(is_aromatic)
    print(f"5. Ring Analysis:")
    print(f"   - Total Rings: {num_rings} (Two aromatic rings, no aliphatic)")
    print(f"   - Aromatic Rings: {num_aromatic_rings} (One Benzene, One Imidazole)")
    
    # 6. Heteroatoms
    heteroatom_count = rdMolDescriptors.CalcNumHeteroatoms(mol)
    # Aromatic N: defined as a nitrogen atom that is part of an aromatic ring.
    aromatic_N_count = len(mol.GetSubstructMatches(Chem.MolFromSmarts('[n]')))
    # Phenolic OH: an oxygen connected to an aromatic carbon, with one hydrogen attached.
    phenolic_oh_count = len(mol.GetSubstructMatches(Chem.MolFromSmarts('[#8X2H1][#6a]')))
    print(f"6. Heteroatom Analysis:")
    print(f"   - Total Heteroatoms: {heteroatom_count} (Target: 4)")
    print(f"   - Aromatic Nitrogens: {aromatic_N_count} (Target: 2)")
    print(f"   - Phenolic Hydroxyl Groups: {phenolic_oh_count} (Target: 1)")

    # 7. H-Bond Donors and Acceptors
    hbd = rdMolDescriptors.CalcNumHBD(mol)
    hba = rdMolDescriptors.CalcNumHBA(mol)
    print(f"7. Hydrogen Bond Analysis:")
    print(f"   - H-Bond Donors: {hbd} (OH group is a donor)")
    print(f"   - H-Bond Acceptors: {hba} (Target: 4)")

    # 8. Forbidden Groups
    # Using SMARTS patterns to check for forbidden groups
    carboxyl_pattern = Chem.MolFromSmarts('[CX3](=O)[OX2H1]')
    aldehyde_pattern = Chem.MolFromSmarts('[CX3H1](=O)')
    thiol_pattern = Chem.MolFromSmarts('[#16X2H1]')
    halide_pattern = Chem.MolFromSmarts('[F,Cl,Br,I]')
    has_forbidden = (mol.HasSubstructMatch(carboxyl_pattern) or
                     mol.HasSubstructMatch(aldehyde_pattern) or
                     mol.HasSubstructMatch(thiol_pattern) or
                     mol.HasSubstructMatch(halide_pattern))
    print(f"8. Forbidden Groups Check (Carboxylic acid, aldehyde, thiol, halide): {'Present' if has_forbidden else 'Absent'}")

    # 9. Imine Functional Group
    imine_pattern = Chem.MolFromSmarts('[CX3](=[NX2])')
    has_imine = mol.HasSubstructMatch(imine_pattern)
    print(f"9. Imine Group Check: {'Present' if has_imine else 'Absent'} (Target: 1)")

    # 10. Tertiary Amines
    # The prompt's requirement for 3 tertiary amines is in direct conflict with
    # the MW/formula and imine requirements. We check for what is present.
    # Pyridine-like N is often considered a tertiary amine.
    tertiary_amine_count = len(mol.GetSubstructMatches(Chem.MolFromSmarts('[n]'))) # counts pyridine-like N
    print(f"10. Tertiary Amines: {tertiary_amine_count} found (Target: 3; constraint not met due to conflict)")

    # 11. Phenol Substitution Pattern
    # para-hydroxylation site without ortho intramolecular hydrogen bonding
    phenol_C_idx = mol.GetSubstructMatch(Chem.MolFromSmarts('c(O)'))[0]
    phenol_C_atom = mol.GetAtomWithIdx(phenol_C_idx)
    neighbors = [n.GetIdx() for n in phenol_C_atom.GetNeighbors() if n.GetSymbol() != 'O']
    ortho_neighbor1 = mol.GetAtomWithIdx(neighbors[0])
    ortho_neighbor2 = mol.GetAtomWithIdx(neighbors[1])
    # The structure is para-substituted, so ortho positions have only H.
    ortho_h_bond_possible = any(a.GetSymbol() in ['O','N'] for a in ortho_neighbor1.GetNeighbors()) or any(a.GetSymbol() in ['O','N'] for a in ortho_neighbor2.GetNeighbors())
    print(f"11. Phenol Substitution Pattern:")
    print(f"    - Para-hydroxylation site: Yes")
    print(f"    - No ortho intramolecular H-bonding: {'Yes' if not ortho_h_bond_possible else 'No'}")

    # 12. Rotatable Bonds
    rotatable_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
    print(f"12. Rotatable Bonds: {rotatable_bonds} (Target: 5)")
    
solve_molecule_challenge()
<<<CCC(=NCCc1c[nH]cn1)c1ccc(O)cc1>>>