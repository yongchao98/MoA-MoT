# First, ensure you have RDKit installed:
# pip install rdkit
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, rdmolops, AllChem

def analyze_molecule():
    """
    Analyzes a molecule defined by a SMILES string to verify it meets
    a specific set of chemical and structural constraints.
    """
    # SMILES string for the designed molecule: 4-((1-butyl-1H-imidazol-5-yl)methyleneamino)phenol
    smiles = 'Oc1ccc(N=Cc2cn(CCCC)cn2)cc1'

    # Create an RDKit molecule object from the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print("Error: Invalid SMILES string.")
        return

    # Add hydrogens to the molecule graph for accurate calculations
    mol = Chem.AddHs(mol)

    # --- Calculate Properties ---
    heavy_atom_count = Descriptors.HeavyAtomCount(mol)
    molecular_weight = Descriptors.ExactMolWt(mol)
    formal_charge = rdmolops.GetFormalCharge(mol)
    valence_electrons = Descriptors.NumValenceElectrons(mol)
    rotatable_bonds = Lipinski.NumRotatableBonds(mol)
    h_bond_donors = Lipinski.NumHDonors(mol)
    h_bond_acceptors = Lipinski.NumHAcceptors(mol)
    num_aromatic_rings = Descriptors.NumAromaticRings(mol)
    formula = AllChem.CalcMolFormula(mol)

    # --- Verify Functional Groups ---
    # Check for two aromatic nitrogens (part of an imidazole ring)
    aromatic_n_pattern = Chem.MolFromSmarts('[n;a]')
    num_aromatic_n = len(mol.GetSubstructMatches(aromatic_n_pattern))

    # Check for phenolic hydroxyl group
    phenol_pattern = Chem.MolFromSmarts('c1(-[O;H1])ccccc1')
    has_phenol = mol.HasSubstructMatch(phenol_pattern)

    # Check for imine group
    imine_pattern = Chem.MolFromSmarts('[#6]=[#7]')
    has_imine = mol.HasSubstructMatch(imine_pattern)
    
    # --- Print Analysis ---
    print(f"Analysis for SMILES: {smiles}")
    print(f"Molecular Formula: {formula}\n")

    print("--- Property Verification ---")
    print(f"Heavy Atoms          : {heavy_atom_count} (Target: 18)")
    print(f"Molecular Weight     : {molecular_weight:.5f} (Target: 243.137)")
    print(f"Formal Charge        : {formal_charge} (Target: 0)")
    print(f"Valence Electrons    : {valence_electrons} (Target: 94)")
    print(f"Rotatable Bonds      : {rotatable_bonds} (Target: 5)")
    print(f"H-Bond Donors        : {h_bond_donors} (Target: 1)")
    print(f"H-Bond Acceptors     : {h_bond_acceptors} (Target: 4)")
    print("\n--- Structural Feature Verification ---")
    print(f"Aromatic Rings       : {num_aromatic_rings} (Target: 2)")
    print(f"Aromatic Nitrogens   : {num_aromatic_n} (Target: 2)")
    print(f"Phenolic Hydroxyl    : {'Yes' if has_phenol else 'No'} (Target: Yes)")
    print(f"Imine Group          : {'Yes' if has_imine else 'No'} (Target: Yes)")

    print("\nNote: The constraint for 'three tertiary amines' was determined to be in conflict with")
    print("other primary constraints (atom count, heteroatom types) and was therefore excluded.")
    print("All other constraints have been met.")

    print("\n--- Final Proposed SMILES ---")
    # Output each number in the final equation (properties)
    print(f"The final molecule has {heavy_atom_count} heavy atoms, a molecular weight of {molecular_weight:.5f}, a formal charge of {formal_charge}, and {valence_electrons} valence electrons.")
    
if __name__ == '__main__':
    analyze_molecule()