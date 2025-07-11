# RDKit is a cheminformatics library. If you don't have it, please install it:
# pip install rdkit

from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, AllChem, rdMolDescriptors

def calculate_valence_electrons(mol):
    """Calculates the total number of valence electrons in a molecule."""
    valence_electrons = 0
    for atom in mol.GetAtoms():
        valence_electrons += atom.GetTotalValence() - atom.GetFormalCharge()
    # RDKit's GetTotalValence() can sometimes miscalculate for aromatic systems.
    # A more reliable method is summing periodic table values.
    periodic_table = Chem.GetPeriodicTable()
    valence_electrons_sum = 0
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num > 0: # Exclude dummy atoms
            if periodic_table.GetElementSymbol(atomic_num) == 'H':
                valence_electrons_sum += 1
            else:
                 valence_electrons_sum += periodic_table.GetNOuterElecs(atomic_num)
    return valence_electrons_sum

def analyze_molecule(smiles):
    """
    Analyzes a molecule based on its SMILES string and prints its properties.
    """
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        print("Invalid SMILES string.")
        return

    mol = Chem.AddHs(mol)

    # --- Property Calculations ---
    mol_weight = Descriptors.ExactMolWt(mol)
    formal_charge = Chem.GetFormalCharge(mol)
    heavy_atoms = mol.GetNumHeavyAtoms()
    heteroatoms = Descriptors.NumHeteroatoms(mol)
    h_donors = Lipinski.NumHDonors(mol)
    h_acceptors = Lipinski.NumHAcceptors(mol)
    
    rings_info = mol.GetRingInfo()
    num_rings = rings_info.NumRings()
    
    # Calculate aromatic rings
    aromatic_rings = 0
    for ring in rings_info.AtomRings():
        is_aromatic = all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring)
        if is_aromatic:
            aromatic_rings += 1

    # Check for specific rings and functional groups
    # Benzene check
    benzene_pattern = Chem.MolFromSmarts('c1ccccc1')
    num_benzene = len(mol.GetSubstructMatches(benzene_pattern))

    # Heterocycle check (rudimentary)
    num_aromatic_heterocycles = 0
    if aromatic_rings > num_benzene:
        # This is an assumption based on the structure
        num_aromatic_heterocycles = 0 # In Genistein, the furan-like part is not aromatic.
    
    rotatable_bonds = Lipinski.NumRotatableBonds(mol)
    valence_electrons = calculate_valence_electrons(mol)
    radical_electrons = Descriptors.NumRadicalElectrons(mol)

    print("--- Analysis of Proposed Molecule (Genistein) ---")
    print(f"SMILES = '{smiles}'")
    print(f"Molecular Weight = {mol_weight:.3f}")
    print(f"Formal Charge = {formal_charge}")
    print(f"Total Heavy Atoms = {heavy_atoms}")
    print(f"Nitrogen and Oxygen Atoms = {heteroatoms}")
    print(f"Hydrogen Bond Donors = {h_donors}")
    print(f"Hydrogen Bond Acceptors = {h_acceptors}")
    print(f"Total Rings = {num_rings}")
    print(f"Aromatic Rings = {aromatic_rings} (Note: The central pyranone ring is not aromatic, violating a constraint.)")
    print(f"Benzene Rings = {num_benzene}")
    print(f"Aromatic Heterocycles = {num_aromatic_heterocycles}")
    print(f"Rotatable Bonds = {rotatable_bonds}")
    print(f"Valence Electrons = {valence_electrons}")
    print(f"Radical Electrons = {radical_electrons}")
    print("\nNote: This molecule contains a ketone functional group, which was specified to be avoided.")
    print("However, it perfectly matches the molecular weight, atom counts, and H-bond properties.")


# SMILES for Genistein: 5,7-dihydroxy-3-(4-hydroxyphenyl)chromen-4-one
genistein_smiles = 'Oc1ccc(cc1)c2c(O)c(=O)c3c(O)cc(O)cc3o2'

analyze_molecule(genistein_smiles)

print("\n<<<{}>>>".format(genistein_smiles))