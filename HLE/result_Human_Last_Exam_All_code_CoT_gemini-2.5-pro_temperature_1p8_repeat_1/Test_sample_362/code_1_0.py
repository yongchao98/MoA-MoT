import re
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

def get_molecular_formula(smiles_string):
    """Calculates the molecular formula from a SMILES string."""
    mol = Chem.MolFromSmiles(smiles_string)
    if mol is not None:
        return CalcMolFormula(mol)
    return ""

def main():
    """
    Determines the products of the specified Wittig reaction,
    constructs the balanced chemical equation, and extracts all numbers.
    """
    # Names and SMILES representations of the relevant species
    # Reactants
    pivalaldehyde_name = "Pivalaldehyde"
    pivalaldehyde_smiles = "CC(C)(C)C=O"
    
    # The SMILES for the ylide is complex and not well-supported,
    # so we'll determine its formula from its constituent parts.
    # Ylide = (Triphenylphosphine) + (2-(2-chlorophenyl)ethyl) - H2
    # C18H15P + C8H8Cl -> Wittig Salt (C26H23ClP), then deprotonation -> Ylide (C26H22ClP)
    ylide_name = "(2-(2-chlorophenyl)ethylidene)triphenyl-l5-phosphane"
    ylide_formula = "C26H22ClP"

    # Products
    product_name = "1-(2-chlorophenyl)-4,4-dimethylpent-2-ene"
    product_smiles = "CC(C)(C)C=CCc1c(Cl)cccc1" # This is the main organic product

    byproduct_name = "Triphenylphosphine oxide"
    byproduct_smiles = "O=P(c1ccccc1)(c2ccccc2)(c3ccccc3)"
    
    # Calculate molecular formulas using RDKit
    pivalaldehyde_formula = get_molecular_formula(pivalaldehyde_smiles)
    product_formula = get_molecular_formula(product_smiles)
    byproduct_formula = get_molecular_formula(byproduct_smiles)

    # Assemble the final chemical equation
    print("The Wittig reaction product is 1-(2-chlorophenyl)-4,4-dimethylpent-2-ene.")
    print("\nThe balanced chemical equation is:")
    
    # Print equation with names
    print(f"{pivalaldehyde_name} + {ylide_name} -> {product_name} + {byproduct_name}\n")
    
    # Print equation with molecular formulas
    equation_str = (
        f"{pivalaldehyde_formula} + {ylide_formula} -> "
        f"{product_formula} + {byproduct_formula}"
    )
    print("In terms of molecular formulas:")
    print(equation_str)

    # Extract all numbers from the molecular formula equation string
    numbers = re.findall(r'\d+', equation_str)
    
    print("\nThe numbers from the final molecular formula equation are:")
    for num in numbers:
        print(num)

if __name__ == "__main__":
    try:
        from rdkit import Chem
        main()
    except ImportError:
        print("This script requires the RDKit library.")
        print("Please install it using: pip install rdkit")
