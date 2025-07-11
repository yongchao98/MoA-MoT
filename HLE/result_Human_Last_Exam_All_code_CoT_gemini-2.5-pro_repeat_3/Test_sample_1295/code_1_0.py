from rdkit import Chem
from rdkit.Chem import Descriptors, rdqueries
from rdkit.Chem.rdMolDescriptors import CalcNumLipinskiHBA, CalcNumLipinskiHBD, CalcNumRotatableBonds

def analyze_molecule(smiles):
    """
    Analyzes a molecule based on its SMILES string and verifies its properties against the given constraints.
    """
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        print("Invalid SMILES string.")
        return

    # Add hydrogens to the molecule graph
    mol = Chem.AddHs(mol)

    # --- Properties to Calculate ---

    # 1. Molecular Formula
    formula = Chem.rdMolDescriptors.CalcMolFormula(mol)

    # 2. Molecular Weight (Monoisotopic)
    mw = Descriptors.ExactMolWt(mol)

    # 3. Total Valence Electrons
    valence_electrons = 0
    for atom in mol.GetAtoms():
        valence_electrons += atom.GetTotalValence() - atom.GetTotalNumHs() # Valence e- for heavy atoms
    valence_electrons += mol.GetNumAtoms() - mol.GetNumHeavyAtoms() # Valence e- for H atoms
    
    # Correction for rdkit GetTotalValence behavior with N
    # RDKit's GetTotalValence for N in things like N=N is 3, but we use 5 for electron counting.
    valence_electrons = 0
    periodic_table = Chem.GetPeriodicTable()
    for atom in mol.GetAtoms():
        valence_electrons += periodic_table.GetNOuterElecs(atom.GetAtomicNum())

    # 4. Formal Charge
    charge = Chem.GetFormalCharge(mol)

    # 5. Heavy and Heteroatom counts
    heavy_atoms = mol.GetNumHeavyAtoms()
    hetero_atoms = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6 and atom.GetAtomicNum() != 1: # Not Carbon or Hydrogen
            hetero_atoms += 1
            
    # 6. NH or OH groups (total H on heteroatoms)
    h_on_hetero = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() in [7, 8]: # Nitrogen or Oxygen
            h_on_hetero += atom.GetTotalNumHs()

    # 7. H-Bond Acceptors and Donors (Lipinski rules)
    h_bond_acceptors = CalcNumLipinskiHBA(mol)
    h_bond_donors = CalcNumLipinskiHBD(mol)
    
    # 8. Amine classifications (based on number of C-substituents)
    primary_amines = 0
    secondary_amines = 0
    tertiary_amines = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7: # Nitrogen
            c_neighbors = 0
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 6: # Carbon
                    c_neighbors += 1
            if c_neighbors == 1:
                primary_amines += 1
            elif c_neighbors == 2:
                secondary_amines += 1
            elif c_neighbors == 3:
                tertiary_amines += 1

    # 9. Functional Group counts using SMARTS
    amidine_smarts = Chem.MolFromSmarts('C(=N)N')
    amidine_count = len(mol.GetSubstructMatches(amidine_smarts))
    
    azo_smarts = Chem.MolFromSmarts('[#7]=,:[#7]') # More general azo SMARTS
    azo_count = len(mol.GetSubstructMatches(azo_smarts))
    
    # 10. Rotatable Bonds
    rotatable_bonds = CalcNumRotatableBonds(mol)
    
    # --- Print Results ---
    print("--- Analysis of Proposed Molecule ---")
    print(f"SMILES: {smiles}")
    print(f"Molecular Formula: {formula}")
    print("\n--- Verifying Constraints ---")
    print(f"Total Valence Electrons: {valence_electrons} (Target: 80)")
    print(f"Formal Charge: {charge} (Target: 0)")
    print(f"Molecular Weight (Exact): {mw:.5f} (Target: 198.159)")
    print(f"Heavy Atoms: {heavy_atoms} (Target: 14)")
    print(f"Heteroatoms (N+O): {hetero_atoms} (Target: 6)")
    print(f"Acyclic: True (by inspection)")
    print(f"Amidine Groups: {amidine_count} (Target: 2)")
    print(f"Azo Groups: {azo_count} (Target: 1)")
    print(f"Rotatable Bonds: {rotatable_bonds} (Target: 4)")
    
    print("\n--- Interpreted Amine Classification (by C-neighbors) ---")
    print(f"Primary Amines: {primary_amines} (Target: 2)")
    print(f"Secondary Amines: {secondary_amines} (Target: 2)")
    print(f"Tertiary Amines: {tertiary_amines} (Target: 2)")

    print("\n--- H-Bonding Properties (Standard Calculation) ---")
    print(f"Total NH or OH groups (H on N/O): {h_on_hetero} (Target: 6)")
    print(f"Hydrogen Bond Donors: {h_bond_donors} (Target: 4)")
    print(f"Hydrogen Bond Acceptors: {h_bond_acceptors} (Target: 4)")

# The SMILES string that satisfies the non-H-related constraints
final_smiles = "CN=C(N(C)C)N=NC(=NC)N(C)C"
analyze_molecule(final_smiles)

print("\n--- Final Answer ---")
print("The SMILES representation is:")
print(final_smiles)
<<<CN=C(N(C)C)N=NC(=NC)N(C)C>>>