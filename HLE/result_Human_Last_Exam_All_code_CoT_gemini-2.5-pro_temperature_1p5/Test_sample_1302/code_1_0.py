import sys
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, Lipinski

def solve_molecule_challenge():
    """
    This function generates and analyzes a molecule based on a specific set of rules.
    Due to contradictory constraints in the prompt (specifically between molecular weight and valence electrons),
    this solution prioritizes the structural rules and the specified molecular weight.
    """
    # SMILES for 4,7-dihydroxy-2-(4-hydroxyphenyl)-5,6-dimethyl-1H-benzo[d]imidazole
    # This structure is chosen as it fits the majority of the complex constraints.
    smiles = 'Cc1c(C)c(O)c2[nH]c(nc2c1O)c3ccc(O)cc3'
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)

    if not mol:
        print("Error: Could not generate molecule from SMILES.")
        return

    # --- Calculations ---
    # Molecular Formula
    formula = rdMolDescriptors.CalcMolFormula(mol)
    
    # Formal Charge
    formal_charge = Chem.GetFormalCharge(mol)

    # Total Molecular Weight (Monoisotopic)
    molecular_weight = Descriptors.MolWt(mol) # Using average mass for comparison
    exact_mw = Descriptors.ExactMolWt(mol)

    # Atom Counts
    heavy_atoms = mol.GetNumHeavyAtoms()
    heteroatoms = Descriptors.NumHeteroatoms(mol)

    # Phenolic Hydroxyl Groups (SMARTS pattern for [OH] attached to an aromatic carbon)
    phenolic_hydroxyl_pattern = Chem.MolFromSmarts('[#8X2H1]-c')
    phenolic_hydroxyls = len(mol.GetSubstructMatches(phenolic_hydroxyl_pattern))

    # Hydrogen Bond Donors and Acceptors
    h_donors = Lipinski.NumHDonors(mol)
    h_acceptors = Lipinski.NumHAcceptors(mol)

    # Ring Information
    ssr = Chem.GetSymmSSSR(mol)
    total_rings = len(ssr)
    aromatic_rings = Lipinski.NumAromaticRings(mol)
    # Benzene rings are aromatic C6 rings
    benzene_pattern = Chem.MolFromSmarts('c1ccccc1')
    benzene_rings = len(mol.GetSubstructMatches(benzene_pattern))
    # Aromatic heterocycles
    aromatic_heterocycles = Lipinski.NumAromaticHeterocycles(mol)

    # Aliphatic/Saturated Rings
    aliphatic_rings = Lipinski.NumAliphaticRings(mol)
    saturated_rings = Lipinski.NumSaturatedRings(mol)

    # Rotatable Bonds
    rotatable_bonds = Descriptors.NumRotatableBonds(mol)

    # Electron Counts
    valence_electrons = Descriptors.NumValenceElectrons(mol)
    # A molecule with 0 formal charge and all paired electrons has 0 radical electrons
    radical_electrons = Descriptors.NumRadicalElectrons(mol)

    # --- Print Results ---
    print(f"Designed Molecule SMILES: {smiles}")
    print("\n--- Properties Verification ---")
    print(f"Formal Charge: {formal_charge}")
    print(f"Molecular Weight (Avg): {molecular_weight:.3f} (Exact: {exact_mw:.3f})")
    print(f"Heavy Atoms: {heavy_atoms}")
    print(f"Heteroatoms (N+O): {heteroatoms}")
    print(f"Phenolic Hydroxyl Groups: {phenolic_hydroxyls}")
    print(f"Hydrogen Bond Donors: {h_donors}")
    print(f"Hydrogen Bond Acceptors: {h_acceptors}")
    print(f"Total Rings: {total_rings}")
    print(f"Aromatic Rings: {aromatic_rings}")
    print(f"  - Benzene Rings: {benzene_rings}")
    print(f"  - Aromatic Heterocycles: {aromatic_heterocycles}")
    print(f"Aliphatic or Saturated Rings: {aliphatic_rings > 0 or saturated_rings > 0}")
    print(f"Rotatable Bonds: {rotatable_bonds}")
    print(f"Valence Electrons: {valence_electrons}")
    print(f"Radical Electrons: {radical_electrons}")
    
    # Print final answer in specified format
    # This final answer is the SMILES string of the designed molecule.
    print(f"\n<<<{smiles}>>>")

# To run the function, we need rdkit installed ('pip install rdkit-pypi').
# If the user doesn't have rdkit, we can't run this code.
# The code is provided for the user to execute in an environment where rdkit is available.
solve_molecule_challenge()
