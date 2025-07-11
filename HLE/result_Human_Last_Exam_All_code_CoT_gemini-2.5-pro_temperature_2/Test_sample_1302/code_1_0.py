import rdkit
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem

def calculate_and_display_molecule_properties():
    """
    This function designs a molecule based on a specific set of rules,
    then calculates and prints its properties.
    """
    # Based on the reasoning, the final molecule is 1-methoxy-4,8-dihydroxydibenzo[b,d]furan.
    # SMILES string for the molecule.
    smiles = "COc1cc2c(cc1O)oc3cc(O)ccc3c2"
    molecule = Chem.MolFromSmiles(smiles)

    # --- Verification of Properties ---
    # Molecular Formula
    formula = Descriptors.rdMolDescriptors.CalcMolFormula(molecule)
    
    # Exact Mass (Monoisotopic)
    exact_mass = Descriptors.ExactMolWt(molecule)
    
    # Valence Electrons
    valence_electrons = sum([atom.GetNumOuterElectrons() for atom in molecule.GetAtoms()]) + Descriptors.calcNumAliphaticCarbocycles(molecule)*0 - Descriptors.calcNumSaturatedCarbocycles(molecule)*0
    
    # Heavy Atoms
    heavy_atoms = Descriptors.HeavyAtomCount(molecule)
    
    # Heteroatoms (N+O)
    num_n = smiles.lower().count('n')
    num_o = smiles.lower().count('o')
    heteroatoms = num_n + num_o
    
    # Rings
    ri = molecule.GetRingInfo()
    num_rings = ri.NumRings()
    num_aromatic_rings = Descriptors.calcNumAromaticRings(molecule)

    # Functional Groups
    phenolic_hydroxyls = Descriptors.NumPhenolicOH(molecule)
    h_bond_donors = Descriptors.NumHDonors(molecule)
    h_bond_acceptors = Descriptors.NumHAcceptors(molecule)
    rotatable_bonds = Descriptors.NumRotatableBonds(molecule)

    print(f"Designed Molecule SMILES: {smiles}")
    print("-" * 30)
    print("Property Verification:")
    print(f"Molecular Formula: {formula}")
    print(f"Formal Charge: {Chem.GetFormalCharge(molecule)}")
    print(f"Valence Electrons: {valence_electrons} (Target: 100)")
    print(f"Heavy Atoms: {heavy_atoms} (Target: 20)")
    print(f"Heteroatoms (N+O): {heteroatoms} (Target: 5)")
    print(f"Total Rings: {num_rings} (Target: 3)")
    print(f"Aromatic Rings: {num_aromatic_rings} (Target: 3)")
    print(f"Phenolic Hydroxyl Groups: {phenolic_hydroxyls} (Target: 3)")
    print(f"Hydrogen Bond Donors: {h_bond_donors} (Target: 3)")
    print(f"Hydrogen Bond Acceptors: {h_bond_acceptors} (Target: 5)")
    print(f"Rotatable Bonds: {rotatable_bonds} (Target: 1)")
    print("-" * 30)

    # Outputting the 'equation' as the exact mass calculation
    print("Final Equation (Exact Mass Calculation):")
    
    num_c = formula.split('C')[1].split('H')[0]
    num_h = formula.split('H')[1].split('O')[0]
    num_o = formula.split('O')[1]

    # Monoisotopic masses
    mass_c = 12.000000
    mass_h = 1.007825
    mass_o = 15.994915

    print(f"({num_c} * C) + ({num_h} * H) + ({num_o} * O) = Exact Mass")
    print(f"({num_c} * {mass_c}) + ({num_h} * {mass_h}) + ({num_o} * {mass_o}) = {exact_mass:.5f} Da (Target: ~270.053 Da)")

# Execute the function
calculate_and_display_molecule_properties()