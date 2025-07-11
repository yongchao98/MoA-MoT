def solve_chemistry_problem():
    """
    Identifies the products of the reaction between styrene and tert-butyl peroxybenzoate.
    """
    print("This reaction is an iron-catalyzed oxo-acyloxylation of styrene.")
    print("It proceeds via a radical mechanism, leading to the formation of two regioisomeric products, A and B.")
    print("-" * 50)

    # Product A: Formed by initial addition of tBuO radical, followed by benzoyloxylation.
    product_a_name = "1-(benzoyloxy)-2-(tert-butoxy)-1-phenylethane"
    product_a_smiles = "c1ccccc1C(OC(=O)c2ccccc2)COC(C)(C)C"

    # Product B: Formed by initial addition of PhCOO radical, followed by tert-butoxylation.
    product_b_name = "2-(benzoyloxy)-1-(tert-butoxy)-1-phenylethane"
    product_b_smiles = "c1ccccc1C(OC(C)(C)C)COC(=O)c2ccccc2"

    print("Product A is:")
    print(f"  Name: {product_a_name}")
    print(f"  SMILES String: {product_a_smiles}")
    print("\nStructure of A: Ph-CH(OCOPh)-CH2-OtBu\n")


    print("Product B is:")
    print(f"  Name: {product_b_name}")
    print(f"  SMILES String: {product_b_smiles}")
    print("\nStructure of B: Ph-CH(OtBu)-CH2-OCOPh\n")
    
    print("-" * 50)
    print("These two compounds are the major products observed as two peaks in the GCMS.")

# Run the function to display the solution.
solve_chemistry_problem()