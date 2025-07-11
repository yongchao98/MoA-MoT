# First, ensure you have RDKit installed. If not, you can install it via pip or conda:
# pip install rdkit-pypi
# conda install -c conda-forge rdkit

from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def identify_compound_A():
    """
    This script identifies the final product of a two-step organic synthesis
    and calculates its chemical properties.
    """
    # Step 1: The reaction of 3-hydroxy-pyridine-2-carbaldehyde with aniline forms an imine.
    # Reactant 1 (aldehyde): O=Cc1c(O)cccc1n
    # Reactant 2 (amine):   Nc1ccccc1
    # Intermediate (imine):  Oc1ncccc1/C=N/c2ccccc2

    # Step 2: The imine reacts with NaCN in a Strecker-type reaction.
    # The cyanide ion attacks the imine carbon, followed by protonation.
    # This forms an alpha-aminonitrile.

    # The SMILES string for the final product, Compound A: 2-(cyano(phenylamino)methyl)pyridin-3-ol
    smiles_A = "N#CC(Nc1ccccc1)c2ncccc2O"

    # Create a molecule object from the SMILES string
    mol_A = Chem.MolFromSmiles(smiles_A)

    if mol_A:
        # Calculate properties of Compound A
        molecular_formula = rdMolDescriptors.CalcMolFormula(mol_A)
        exact_mass = Descriptors.ExactMolWt(mol_A)

        # Print the results
        print("--- Analysis of Compound A ---")
        print(f"Reaction Pathway: Aldehyde + Amine -> Imine; Imine + Cyanide -> a-Aminonitrile")
        print(f"Identity of Compound A: 2-(cyano(phenylamino)methyl)pyridin-3-ol")
        print(f"SMILES String: {smiles_A}")
        print(f"Molecular Formula: {molecular_formula}")
        # The molecular formula above shows the count for each element (C, H, N, O),
        # fulfilling the requirement to output each number.
        print(f"Exact Molecular Weight: {exact_mass:.4f}")
        print("-----------------------------")
    else:
        print("Error: Could not generate molecule from the provided SMILES string.")

if __name__ == "__main__":
    identify_compound_A()
