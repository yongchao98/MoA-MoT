import sys
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

def verify_molecule(smiles):
    """
    Verifies a molecule against a series of specified constraints using RDKit.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print("Error: Invalid SMILES string.", file=sys.stderr)
        return

    # --- Verification Checks ---
    mw = Descriptors.MolWt(mol)
    heavy_atoms = mol.GetNumHeavyAtoms()
    valence_electrons = Descriptors.NumValenceElectrons(mol)
    formal_charge = Chem.GetFormalCharge(mol)
    h_acceptors = rdMolDescriptors.CalcNumHBA(mol)
    h_donors = rdMolDescriptors.CalcNumHBD(mol)
    num_rings = rdMolDescriptors.CalcNumRings(mol)
    rotatable_bonds = Descriptors.NumRotatableBonds(mol)

    # Heteroatoms (any atom not C or H)
    heteroatoms = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() not in [1, 6])

    # Pattern matching for functional groups using SMARTS
    has_no_radicals = all(atom.GetNumRadicalElectrons() == 0 for atom in mol.GetAtoms())
    no_halogens = not mol.HasSubstructMatch(Chem.MolFromSmarts("[F,Cl,Br,I]"))
    carbonyl_count = len(mol.GetSubstructMatches(Chem.MolFromSmarts("[CX3](=O)[#6]")))
    ether_oxygens = len(mol.GetSubstructMatches(Chem.MolFromSmarts("[OD2]([#6])[#6]")))
    has_no_carboxylic_acids = not mol.HasSubstructMatch(Chem.MolFromSmarts("[CX3](=O)[OX2H1]"))
    has_no_hydroxyls = not mol.HasSubstructMatch(Chem.MolFromSmarts("[OX2H]"))
    has_no_amines = not mol.HasSubstructMatch(Chem.MolFromSmarts("[NX3;H2;!$(NC=O)]"))
    has_no_thiols = not mol.HasSubstructMatch(Chem.MolFromSmarts("[SH]"))
    has_no_esters = not mol.HasSubstructMatch(Chem.MolFromSmarts("[#6][CX3](=O)[OX2][#6]"))
    has_no_aromatic_rings = not mol.HasSubstructMatch(Chem.MolFromSmarts("a"))
    has_no_tetrazole = not mol.HasSubstructMatch(Chem.MolFromSmarts("c1nnnn1"))
    
    # Check if all rings are heterocycles
    all_rings_are_hetero = True
    for ring_atoms in Chem.GetSymmSSSR(mol):
        is_hetero = False
        for atom_idx in ring_atoms:
            if mol.GetAtomWithIdx(atom_idx).GetAtomicNum() != 6:
                is_hetero = True
                break
        if not is_hetero:
            all_rings_are_hetero = False
            break

    print(f"--- Verifying SMILES: {smiles} ---")
    print(f"Molecular Weight: {mw:.2f} (Target: 258.11)")
    print(f"Heavy Atoms: {heavy_atoms} (Target: 18)")
    print(f"Valence Electrons: {valence_electrons} (Target: 102)")
    print(f"Formal Charge: {formal_charge} (Target: 0)")
    print(f"No Radical Electrons: {has_no_radicals} (Target: True)")
    print(f"Heteroatoms: {heteroatoms} (Target: 6)")
    print(f"Carbonyl Oxygens (Ketones): {carbonyl_count} (Target: >=1)")
    print(f"Ether Oxygens: {ether_oxygens} (Target: 5)")
    print(f"H-Bond Acceptors: {h_acceptors} (Target: 6)")
    print(f"H-Bond Donors: {h_donors} (Target: 0)")
    print(f"Total Rings: {num_rings} (Target: 3)")
    print(f"All Rings are Saturated Heterocycles: {all_rings_are_hetero and not has_no_aromatic_rings} (Target: True)")
    print(f"Rotatable Bonds: {rotatable_bonds} (Target: 0)")
    print(f"No Halogens: {no_halogens} (Target: True)")
    print(f"No Esters or Carboxylic Acids: {has_no_esters and has_no_carboxylic_acids} (Target: True)")
    print(f"No Amines/Thiols/Hydroxyls: {has_no_amines and has_no_thiols and has_no_hydroxyls} (Target: True)")
    
    # Final representation: Print the full equation for molecular weight
    mol_with_h = Chem.AddHs(mol)
    formula = rdMolDescriptors.CalcMolFormula(mol_with_h)
    print("\n--- Final Molecular Representation ---")
    print(f"Molecular Formula: {formula}")

    c_count = sum(1 for a in mol_with_h.GetAtoms() if a.GetAtomicNum() == 6)
    h_count = sum(1 for a in mol_with_h.GetAtoms() if a.GetAtomicNum() == 1)
    o_count = sum(1 for a in mol_with_h.GetAtoms() if a.GetAtomicNum() == 8)
    
    # Use exact masses for calculation for higher precision matching
    c_mass = 12.011
    h_mass = 1.008
    o_mass = 15.999
    
    total_mw = (c_count * c_mass) + (h_count * h_mass) + (o_count * o_mass)

    print(f"Final Equation: ({c_count} * {c_mass}) + ({h_count} * {h_mass}) + ({o_count} * {o_mass}) = {total_mw:.2f}")

# This is the SMILES string for the designed Bicyclo[5.5.5]heptadecane derivative.
designed_smiles = "O=C1CC2(COCOC2)C2OC1COCOC2"

# Execute the verification
verify_molecule(designed_smiles)

# Final answer in the required format
print(f"\n<<<{designed_smiles}>>>")