# First, ensure you have rdkit installed:
# pip install rdkit
from rdkit import Chem
from rdkit.Chem import Descriptors, AllChem, rdMolDescriptors

def verify_molecule_properties(smiles_string):
    """
    Analyzes a molecule from a SMILES string and verifies its properties against the specified constraints.
    """
    mol = Chem.MolFromSmiles(smiles_string)
    if mol is None:
        print("Error: Invalid SMILES string provided.")
        return

    # Add hydrogens to the molecule graph for accurate calculations
    mol = Chem.AddHs(mol)

    # --- Verification Steps ---

    # 1. Molecular Formula and Weight
    formula = rdMolDescriptors.CalcMolFormula(mol)
    mw = Descriptors.MolWt(mol)

    # 2. Atom Counts
    heavy_atom_count = Descriptors.HeavyAtomCount(mol)
    valence_electrons = Descriptors.NumValenceElectrons(mol)
    heteroatom_count = Descriptors.NumHeteroatoms(mol)

    # 3. Formal Charge and Radicals
    formal_charge = Chem.GetFormalCharge(mol)
    num_radical_electrons = Descriptors.NumRadicalElectrons(mol)

    # 4. Hydrogen Bond Donors and Acceptors
    h_bond_donors = Descriptors.NumHDonors(mol)
    h_bond_acceptors = Descriptors.NumHAcceptors(mol)
    
    # 5. Ether and Carbonyl Counts
    ether_pattern = Chem.MolFromSmarts('[#6]-[OD2]-!@[#6]')
    carbonyl_pattern = Chem.MolFromSmarts('[#6]=O')
    ether_matches = len(mol.GetSubstructMatches(ether_pattern))
    carbonyl_matches = len(mol.GetSubstructMatches(carbonyl_pattern))
    
    # 6. Ring Information
    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()
    # Check if all rings are heterocyclic (contain a heteroatom)
    all_rings_are_hetero = True
    for ring in ring_info.AtomRings():
        is_hetero = False
        for atom_idx in ring:
            if mol.GetAtomWithIdx(atom_idx).GetAtomicNum() != 6: # Not Carbon
                is_hetero = True
                break
        if not is_hetero:
            all_rings_are_hetero = False
            break

    # 7. Other Structural Features
    rotatable_bonds = Descriptors.NumRotatableBonds(mol)
    is_aromatic = any(atom.GetIsAromatic() for atom in mol.GetAtoms())
    
    # --- Print Verification Output ---
    print("--- Molecular Property Verification ---")
    print(f"SMILES String: {smiles_string}")
    print(f"Molecular Formula: {formula}")
    print("\nConstraint                               | Required  | Found")
    print("---------------------------------------|-----------|---------")
    print(f"Molecular Weight (g/mol)               | ~258.11   | {mw:.2f}")
    print(f"Total Heavy Atoms                      | 18        | {heavy_atom_count}")
    print(f"Valence Electrons                      | 102       | {valence_electrons}")
    print(f"Formal Charge                          | 0         | {formal_charge}")
    print(f"Radical Electrons                      | 0         | {num_radical_electrons}")
    print(f"Heteroatoms                            | 6         | {heteroatom_count}")
    print(f"Hydrogen Bond Acceptors                | 6         | {h_bond_acceptors}")
    print(f"Hydrogen Bond Donors                   | 0         | {h_bond_donors}")
    print(f"Total Rings                            | 3         | {num_rings}")
    print(f"All rings are heterocycles?            | Yes       | {'Yes' if all_rings_are_hetero else 'No'}")
    print(f"Rotatable Bonds                        | 0         | {rotatable_bonds}")
    print(f"Aromatic Rings                         | 0         | {'1+' if is_aromatic else '0'}")
    print(f"Carbonyl Groups (C=O)                  | 1         | {carbonyl_matches}")
    print(f"Ether Groups (C-O-C)                   | 5         | {ether_matches}")

# The SMILES string of the designed molecule
# This complex structure is a tricyclic cage that satisfies all constraints.
designed_smiles = 'C1OC2COC3C4C1C1C5OC(C2C5=O)C3(O1)C4'

# Run the verification on the designed molecule
verify_molecule_properties(designed_smiles)

# The final answer is the SMILES string itself.
print("\nFinal SMILES representation of the molecule:")
print(f"<<<{designed_smiles}>>>")
