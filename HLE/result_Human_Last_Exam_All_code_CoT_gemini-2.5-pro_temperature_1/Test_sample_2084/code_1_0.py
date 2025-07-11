# First, ensure you have RDKit installed:
# pip install rdkit

from rdkit import Chem
from rdkit.Chem import Descriptors

def solve_reaction_product_mass():
    """
    This script identifies the products of a chemical reaction and determines which has the higher molar mass.
    """
    print("Step 1: Analyzing the reaction.")
    print("The reaction involves a complex molecule in an acidic aqueous solution (TFA in water).")
    print("The provided SMILES 'CC12COC(OC1)(OC2)C1=CC=CC=C1' is invalid, but its fragments suggest a spiroketal structure.")
    print("Spiroketals undergo acid-catalyzed hydrolysis to form a ketone and a diol.\n")

    print("Step 2: Deducing the products based on the fragments.")
    # The 'CC' (methyl) and 'C1=CC=CC=C1' (phenyl) fragments point to Acetophenone.
    # The 'COC' (ether) fragment points to a diol like Diethylene Glycol.
    
    # Product 1: Acetophenone
    ketone_name = "Acetophenone"
    ketone_smiles = "CC(=O)c1ccccc1"
    
    # Product 2: Diethylene Glycol
    diol_name = "Diethylene Glycol"
    diol_smiles = "OOCCOCCO"

    print(f"Product 1 (Ketone) is likely {ketone_name}.")
    print(f"Product 2 (Diol) is likely {diol_name}.\n")

    print("Step 3: Calculating molar masses.")
    # Create RDKit molecule objects
    ketone_mol = Chem.MolFromSmiles(ketone_smiles)
    diol_mol = Chem.MolFromSmiles(diol_smiles)

    # Calculate exact molecular weight
    ketone_mw = Descriptors.ExactMolWt(ketone_mol)
    diol_mw = Descriptors.ExactMolWt(diol_mol)

    print(f"Details of Product 1 ({ketone_name}):")
    print(f"  SMILES: {ketone_smiles}")
    print(f"  Molar Mass: {ketone_mw:.4f} g/mol")

    print(f"Details of Product 2 ({diol_name}):")
    print(f"  SMILES: {diol_smiles}")
    print(f"  Molar Mass: {diol_mw:.4f} g/mol\n")

    print("Step 4: Comparing the molar masses.")
    if ketone_mw > diol_mw:
        higher_mass_product_name = ketone_name
        higher_mass_product_smiles = ketone_smiles
    else:
        higher_mass_product_name = diol_name
        higher_mass_product_smiles = diol_smiles

    print(f"The product with the higher molar mass is {higher_mass_product_name}.")
    print(f"Final Answer (SMILES string): {higher_mass_product_smiles}")

# Execute the function
solve_reaction_product_mass()