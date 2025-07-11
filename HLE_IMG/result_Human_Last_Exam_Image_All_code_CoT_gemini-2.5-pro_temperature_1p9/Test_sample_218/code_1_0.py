# To run this code, you need to install the rdkit library:
# pip install rdkit

from rdkit import Chem
from rdkit.Chem import Descriptors

def get_molecule_info(name, smiles_string):
    """Creates a molecule from SMILES and returns its formula and weight."""
    molecule = Chem.MolFromSmiles(smiles_string)
    if molecule:
        formula = Descriptors.rdMolDescriptors.CalcMolFormula(molecule)
        mol_weight = Descriptors.ExactMolWt(molecule)
        print(f"Molecule: {name}")
        print(f"  SMILES: {smiles_string}")
        print(f"  Formula: {formula}")
        print(f"  Molecular Weight: {mol_weight:.4f}\n")
    else:
        print(f"Could not parse SMILES for {name}")

# --- Analysis ---
print("Analyzing the deoxygenation of Geraniol to Compound A:\n")

# Geraniol ((2E)-3,7-dimethylocta-2,6-dien-1-ol)
# Note: RDKit SMILES do not always explicitly show stereochemistry in this simple format.
geraniol_smiles = "CC(=CCCC(=CCO)C)C"
get_molecule_info("Geraniol (Starting Material)", geraniol_smiles)

# The product, Compound A, is formed by replacing the -OH group with an -H.
# Structure: (CH3)2C=CH-CH2-CH2-C(CH3)=CH-CH3
# IUPAC Name: 2,6-dimethylocta-2,6-diene
product_A_smiles = "CC(=CCCC(=CC)C)C"
get_molecule_info("Compound A (Final Product)", product_A_smiles)

print("The transformation from C10H18O (Geraniol) to C10H18 (Compound A) confirms the deoxygenation reaction.")
print("The final product, Compound A, is 2,6-dimethylocta-2,6-diene.")
