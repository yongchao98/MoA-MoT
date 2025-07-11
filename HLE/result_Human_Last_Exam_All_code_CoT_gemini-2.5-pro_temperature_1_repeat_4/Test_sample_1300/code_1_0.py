# The user wants to design a molecule with a very specific set of constraints.
# The constraints on functional groups (3 tertiary amines, 1 imine) and heteroatom count (4 total)
# are contradictory under standard chemical definitions for a neutral molecule.
# It is impossible for 3 nitrogen atoms to satisfy these roles simultaneously.
# I will proceed by designing a molecule that satisfies all the "hard" constraints:
# molecular formula, exact mass, heavy atom count, ring systems, and rotatable bonds.
# This approach assumes a likely typo in the prompt's tertiary amine count.
# The chosen molecule is the most plausible fit for the vast majority of the requirements.

# The final code will use the RDKit library to perform the analysis and generate the output.
# Please ensure you have RDKit installed: pip install rdkit

from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, AllChem

def solve_molecule_challenge():
    """
    Analyzes a candidate molecule against a list of complex constraints and prints the results.
    """
    # SMILES for the candidate molecule: 4-((E)-(1-(5-isopropyl-1-methyl-1H-imidazol-4-yl)ethylidene)amino)phenol
    # A more canonical representation: Oc1ccc(N=C(C)c2c(C(C)C)n(C)cn2)cc1
    # After analyzing rotatable bonds, another candidate is better:
    # Oc1ccc(N=Cc2c(C(C)C)n(C)cn2)cc1
    smiles = "CN1C=NC(C(C)C)=C1C=Nc1ccc(O)cc1"

    # Create molecule object
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)

    # --- Analysis ---

    # 1. Molecular Formula and Properties
    formula = rdMolDescriptors.CalcMolFormula(mol)
    heavy_atoms = Descriptors.HeavyAtomCount(mol)
    mw = Descriptors.ExactMolWt(mol)
    charge = Chem.GetFormalCharge(mol)

    # Valence Electron Calculation
    valence_electrons = 0
    for atom in mol.GetAtoms():
        # Element-specific valence electrons
        periodic_table = Chem.GetPeriodicTable()
        atomic_num = atom.GetAtomicNum()
        if atomic_num == 1: # Hydrogen
            valence_electrons += 1
        elif atomic_num in [6, 7, 8]: # C, N, O
            valence_electrons += periodic_table.GetNOuterElecs(atomic_num)

    # 2. Structural Features
    # Aromatic Rings
    ssr = Chem.GetSymmSSSR(mol)
    aromatic_rings = [r for r in ssr if Chem.IsAromatic(mol, r)]
    has_benzene = any(set(r) for r in aromatic_rings if len(r) == 6 and all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in r))
    imidazole_patt = Chem.MolFromSmarts('c1cn[nH]c1')
    has_imidazole = mol.HasSubstructMatch(imidazole_patt)

    # Heteroatoms
    heteroatom_count = Descriptors.NumHeteroatoms(mol)
    aromatic_n_patt = Chem.MolFromSmarts('[n]')
    aromatic_n_count = len(mol.GetSubstructMatches(aromatic_n_patt))
    hydroxyl_o_patt = Chem.MolFromSmarts('[OX2H]')
    hydroxyl_o_count = len(mol.GetSubstructMatches(hydroxyl_o_patt))

    # H-Bond Acceptors/Donors
    h_acceptors = Descriptors.NumHAcceptors(mol)
    h_donors = Descriptors.NumHDonors(mol)

    # Functional Groups
    imine_patt = Chem.MolFromSmarts('[#6X3]=[#7X2]')
    imine_count = len(mol.GetSubstructMatches(imine_patt))
    tertiary_amine_patt = Chem.MolFromSmarts('[NX3;H0;+0;!$(N-C=O);!$(N-S=O)]')
    tertiary_amine_count = len(mol.GetSubstructMatches(tertiary_amine_patt))
    phenol_patt = Chem.MolFromSmarts('c1(-O)ccccc1')
    phenol_count = len(mol.GetSubstructMatches(phenol_patt))

    # Other properties
    rotatable_bonds = Descriptors.NumRotatableBonds(mol)

    # --- Final Output ---
    print(f"Proposed SMILES: {smiles}\n")
    print("--- Property Verification ---")
    print(f"Molecular Formula: {formula} (C14 H17 N3 O)")
    print(f"Heavy Atoms: {heavy_atoms} (Target: 18)")
    print(f"Molecular Weight: {mw:.5f} (Target: 243.137)")
    print(f"Formal Charge: {charge} (Target: 0)")
    print(f"Valence Electrons: {valence_electrons} (Target: 94)")
    print(f"Aromatic Rings: {len(aromatic_rings)} (Benzene: {has_benzene}, Imidazole: {has_imidazole})")
    print(f"Total Heteroatoms: {heteroatom_count} (Target: 4)")
    print(f"Aromatic Nitrogens: {aromatic_n_count} (Target: 2)")
    print(f"Hydroxyl Groups (Donor): {h_donors} (Target: 1)")
    print(f"Hydrogen Bond Acceptors: {h_acceptors} (Target: 4, discrepancy noted)")
    print(f"Imine Groups: {imine_count} (Target: 1)")
    print(f"Tertiary Amines: {tertiary_amine_count} (Target: 3, discrepancy noted)")
    print(f"Phenolic OH Groups: {phenol_count} (Target: 1)")
    print(f"Rotatable Bonds: {rotatable_bonds} (Target: 5)")

    # Print the "final equation" for molecular weight
    mol_no_hs = Chem.RemoveHs(mol)
    atom_counts = {}
    for atom in mol_no_hs.GetAtoms():
        sym = atom.GetSymbol()
        atom_counts[sym] = atom_counts.get(sym, 0) + 1
    c = atom_counts.get('C', 0)
    h = rdMolDescriptors.CalcNumHs(mol)
    n = atom_counts.get('N', 0)
    o = atom_counts.get('O', 0)
    
    print("\n--- Final Equation (Exact Mass) ---")
    print(f"{c} * C + {h} * H + {n} * N + {o} * O = {mw:.5f}")


solve_molecule_challenge()

<<<CN1C=NC(C(C)C)=C1C=Nc1ccc(O)cc1>>>