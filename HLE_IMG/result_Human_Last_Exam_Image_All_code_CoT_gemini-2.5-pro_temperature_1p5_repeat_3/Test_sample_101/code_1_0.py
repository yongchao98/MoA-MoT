# To run this code, you may need to install the RDKit library:
# pip install rdkit

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def get_molecular_formula(smiles_string):
    """A helper function to calculate the molecular formula from a SMILES string."""
    molecule = Chem.MolFromSmiles(smiles_string)
    if molecule:
        return rdMolDescriptors.CalcMolFormula(molecule)
    return "Invalid SMILES"

def identify_reaction_product():
    """
    Identifies Compound A from the described two-step reaction
    and presents the balanced chemical equation for the overall transformation.
    """
    # Define the SMILES strings for the molecules involved in the overall reaction
    reactant_aldehyde_smiles = "O=Cc1c(O)cccn1"  # 3-hydroxy-pyridine-2-carbaldehyde
    reactant_amine_smiles = "c1ccccc1N"          # Aniline
    product_A_smiles = "N#CC(Nc1ccccc1)c2c(O)cccn2" # (phenylamino)(3-hydroxypyridin-2-yl)acetonitrile
    byproduct_water_smiles = "O"                 # Water

    # Get the molecular formulas for each molecule
    formula_aldehyde = get_molecular_formula(reactant_aldehyde_smiles)
    formula_amine = get_molecular_formula(reactant_amine_smiles)
    # The source of the cyanide group is NaCN. The overall reaction can be conceptually
    # written with HCN, as a proton is also consumed.
    formula_cyanide_source = "CHN"
    formula_A = get_molecular_formula(product_A_smiles)
    formula_water = get_molecular_formula(byproduct_water_smiles)

    print("--- Reaction Summary ---")
    print("The reaction is a Strecker synthesis, which transforms an aldehyde, an amine, and a cyanide source into an alpha-aminonitrile.")

    print("\n--- Overall Balanced Chemical Equation ---")
    # The stoichiometry of the overall reaction is 1:1:1 -> 1:1.
    # The prompt requires printing each number in the equation.
    print(f"1 {formula_aldehyde} (Aldehyde) + 1 {formula_amine} (Aniline) + 1 {formula_cyanide_source} (from NaCN) -> 1 {formula_A} (Compound A) + 1 {formula_water} (Water)")

    print("\n--- Compound A Details ---")
    molecule_A = Chem.MolFromSmiles(product_A_smiles)
    if molecule_A:
        molecular_weight_A = rdMolDescriptors.CalcExactMolWt(molecule_A)
        print(f"Chemical Name: (phenylamino)(3-hydroxypyridin-2-yl)acetonitrile")
        print(f"Molecular Formula: {formula_A}")
        print(f"SMILES: {product_A_smiles}")
        print(f"Exact Molecular Weight: {molecular_weight_A:.4f} g/mol")
    else:
        print("Could not generate details for Compound A.")

if __name__ == "__main__":
    identify_reaction_product()