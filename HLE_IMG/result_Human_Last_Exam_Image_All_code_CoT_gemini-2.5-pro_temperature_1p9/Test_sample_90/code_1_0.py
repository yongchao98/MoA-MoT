def get_product_structures():
    """
    This function describes the chemical structures of products A, B, and C.
    """
    
    product_A = {
        "name": "Product A",
        "formula": "C14H20N2O3",
        "structure_description": (
            "A substituted pyrrolizidine. It is the result of the [3+2] cycloaddition "
            "between the ylide (from decarboxylation of the starting material) and methyl propiolate, "
            "followed by N-acylation of the pendant dihydropyrrole ring. It contains a methyl ester group and an N-acetyl group."
        ),
        "smiles_representation": "COC(=O)C1=C[C@H]2[C@@H](C1)N(CCC2)C3=N(C(=O)C)C=CC3"
    }

    product_B = {
        "name": "Product B",
        "formula": "C12H14N2O3",
        "structure_description": (
            "A heavily oxidized and dehydrogenated derivative of the initial cycloadduct. It features a conjugated "
            "pyrrolizidinone core (a bicyclic lactam) and a pyrrole substituent, which accounts for the downfield NMR signals. "
            "The structure results from oxidation of the saturated pyrrolizidine core and dehydrogenation of the heterocyclic rings."
        ),
        "smiles_representation": "COC(=O)C1=C[C@@H]2N(C(=O)CC2)c3nccc13" # Example isomer, exact structure can vary
    }
    
    product_C = {
        "name": "Product C",
        "formula": "C11H16N2O3",
        "structure_description": (
            "A direct acylation product of the starting material, without undergoing cycloaddition. An acetyl group "
            "has been added to the C3 position of the dihydropyrrole ring. The core structure is "
            "1-(3-acetyl-4,5-dihydro-3H-pyrrol-2-yl)pyrrolidine-2-carboxylic acid."
        ),
        "smiles_representation": "O=C(O)[C@H]1N(CCC1)C2=NC(C(=O)C)CC2"
    }
    
    products = [product_A, product_B, product_C]
    
    for product in products:
        print(f"--- {product['name']} ---")
        print(f"Molecular Formula: {product['formula']}")
        print("\nProposed Structure:")
        print(product['structure_description'])
        # The SMILES strings are complex and represent one possible isomer; they are provided for illustrative purposes.
        # print(f"Illustrative SMILES: {product['smiles_representation']}")
        print("-" * (len(product['name']) + 8) + "\n")

# Execute the function to print the solution
get_product_structures()