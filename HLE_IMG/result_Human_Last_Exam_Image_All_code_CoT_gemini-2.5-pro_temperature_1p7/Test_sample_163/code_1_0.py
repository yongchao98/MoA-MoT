import sys

# Suppress warnings if any library gives them, for cleaner output.
# (Not strictly necessary here but good practice)
if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")

def solve_reaction():
    """
    Identifies the products of the given chemical reaction and prints the results.
    """
    # Molecular formulas of reactants
    styrene_formula = "C8H8"
    peroxide_formula = "C11H14O3"
    
    # Molecular formula of the products (they are isomers)
    product_formula = "C19H22O3"

    print("The overall reaction is an addition, with Fe(OTf)3 as a catalyst.")
    print("The balanced chemical equation is:")
    # Print the equation with stoichiometric numbers (which are all 1)
    print(f"1 {styrene_formula} + 1 {peroxide_formula} --> 1 {product_formula}")
    print("-" * 40)
    print("The reaction yields two major products, A and B, which are regioisomers.\n")
    
    # Information for Product A
    product_A = {
        "name": "2-(tert-butoxy)-1-phenylethyl benzoate",
        "structure_sketch": "Ph-CH(OCOPh)-CH2-OtBu",
        "smiles": "CC(C)(C)OC[CH](OC(=O)c1ccccc1)c2ccccc2"
    }

    # Information for Product B
    product_B = {
        "name": "2-(benzoyloxy)-1-(tert-butoxy)-1-phenylethane",
        "structure_sketch": "Ph-CH(OtBu)-CH2-OCOPh",
        "smiles": "CC(C)(C)OC(c1ccccc1)COC(=O)c2ccccc2"
    }

    # Print product details
    print("Product A:")
    print(f"  Name: {product_A['name']}")
    print(f"  Structure: {product_A['structure_sketch']}")
    print(f"  SMILES: {product_A['smiles']}\n")

    print("Product B:")
    print(f"  Name: {product_B['name']}")
    print(f"  Structure: {product_B['structure_sketch']}")
    print(f"  SMILES: {product_B['smiles']}")

# Execute the function to print the solution
solve_reaction()