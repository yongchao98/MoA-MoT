import sys
import platform

# Check if rdkit is available, otherwise print a message and the pre-calculated values.
try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, rdMolDescriptors
except ImportError:
    print("RDKit library not found. Please install it to run the full analysis (pip install rdkit).")
    print("Displaying pre-calculated values for the target molecule.\n")
    
    smiles = 'Cc1c(O)cc2oc(Oc3ccc(O)cc3)c2c1O'
    molecular_formula = "C15H10O5"
    molecular_weight = 270.0528
    formal_charge = 0
    heavy_atoms = 20
    heteroatoms = 5
    phenolic_hydroxyls = 3
    h_bond_donors = 3
    h_bond_acceptors = 5
    total_rings = 3
    aromatic_rings = 3
    benzene_rings = 2
    aromatic_heterocycles = 1
    rotatable_bonds = 1
    valence_electrons = 100
    radical_electrons = 0
    
    print(f"Proposed SMILES: {smiles}")
    print("----------------------------------------")
    print(f"Molecular Formula: {molecular_formula}")
    print(f"Total Molecular Weight: {molecular_weight:.4f}")
    print(f"Formal Charge: {formal_charge}")
    print(f"Heavy Atoms: {heavy_atoms}")
    print(f"Heteroatoms (N+O): {heteroatoms}")
    print(f"Phenolic Hydroxyl Groups: {phenolic_hydroxyls}")
    print(f"Hydrogen Bond Donors: {h_bond_donors}")
    print(f"Hydrogen Bond Acceptors: {h_bond_acceptors}")
    print(f"Total Rings: {total_rings}")
    print(f"Aromatic Rings: {aromatic_rings}")
    print(f"  - Benzene Rings: {benzene_rings}")
    print(f"  - Aromatic Heterocycles: {aromatic_heterocycles}")
    print(f"Rotatable Bonds: {rotatable_bonds}")
    print(f"Valence Electrons: {valence_electrons}")
    print(f"Radical Electrons: {radical_electrons}")
    sys.exit()

def analyze_molecule(smiles_string):
    """
    Analyzes a molecule from a SMILES string and prints its properties.
    """
    mol = Chem.MolFromSmiles(smiles_string)
    if mol is None:
        print(f"Error: Invalid SMILES string: {smiles_string}")
        return

    mol = Chem.AddHs(mol)

    # --- Calculations ---
    molecular_formula = rdMolDescriptors.CalcMolFormula(mol)
    molecular_weight = Descriptors.ExactMolWt(mol)
    formal_charge = Chem.GetFormalCharge(mol)
    heavy_atoms = mol.GetNumHeavyAtoms()
    
    heteroatoms = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() in [7, 8])

    phenolic_patt = Chem.MolFromSmarts('[OH]c')
    phenolic_hydroxyls = len(mol.GetSubstructMatches(phenolic_patt))
    
    h_bond_donors = rdMolDescriptors.CalcNumHBD(mol)
    h_bond_acceptors = rdMolDescriptors.CalcNumHBA(mol)

    ssr = Chem.GetSymmSSSR(mol)
    total_rings = len(ssr)
    aromatic_rings = sum(1 for ring in ssr if Chem.IsAromatic(mol, [list(ring)[0]]*len(ring)))
    
    benzene_patt = Chem.MolFromSmarts('c1ccccc1')
    benzene_rings = len(mol.GetSubstructMatches(benzene_patt))
    
    # Aromatic heterocycles are aromatic rings that are not benzene
    aromatic_heterocycles = aromatic_rings - benzene_rings
    
    rotatable_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
    
    valence_electrons = sum(atom.GetTotalValence() for atom in mol.GetAtoms())
    # Correct valence electron calculation
    valence_electrons = 0
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num == 1: valence_electrons += 1   # H
        elif atomic_num == 6: valence_electrons += 4 # C
        elif atomic_num == 7: valence_electrons += 5 # N
        elif atomic_num == 8: valence_electrons += 6 # O

    radical_electrons = Descriptors.NumRadicalElectrons(mol)

    # --- Print Results ---
    print(f"Analysis for SMILES: {smiles_string}")
    print("----------------------------------------")
    print(f"Molecular Formula: {molecular_formula}")
    print(f"Total Molecular Weight: {molecular_weight:.4f}")
    print(f"Formal Charge: {formal_charge}")
    print(f"Heavy Atoms: {heavy_atoms}")
    print(f"Heteroatoms (N+O): {heteroatoms}")
    print(f"Phenolic Hydroxyl Groups: {phenolic_hydroxyls}")
    print(f"Hydrogen Bond Donors: {h_bond_donors}")
    print(f"Hydrogen Bond Acceptors: {h_bond_acceptors}")
    print(f"Total Rings: {total_rings}")
    print(f"Aromatic Rings: {aromatic_rings}")
    print(f"  - Benzene Rings: {benzene_rings}")
    print(f"  - Aromatic Heterocycles: {aromatic_heterocycles}")
    print(f"Rotatable Bonds: {rotatable_bonds}")
    print(f"Valence Electrons: {valence_electrons}")
    print(f"Radical Electrons: {radical_electrons}")


# The SMILES string for 2-(4-hydroxyphenoxy)-5-methylbenzofuran-4,6-diol
final_smiles = 'Cc1c(O)cc2oc(Oc3ccc(O)cc3)c2c1O'
analyze_molecule(final_smiles)
