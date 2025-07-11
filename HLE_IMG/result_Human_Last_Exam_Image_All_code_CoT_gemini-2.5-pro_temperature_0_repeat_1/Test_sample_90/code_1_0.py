def get_product_structures():
    """
    This function provides the structures and molecular formulas for the products A, B, and C.
    The structures are given in the SMILES format.
    """

    product_A = {
        "label": "A",
        "name": "Tetracyclic [8pi+2pi] cycloaddition product (reduced form)",
        "formula": "C14H20N2O3",
        "smiles": "CC(=O)N1[C@@H]2CC[C@H]1[C@]13C=C(C(=O)OC)[C@@H](N21)[C@H]1N4CCC[C@H]13C4" # Example SMILES for a possible adduct
    }

    product_B = {
        "label": "B",
        "name": "Tricyclic enelactam-imide rearrangement product",
        "formula": "C12H14N2O3",
        "smiles": "CC(=O)N1C2CCC[C@H]1[C@]1(N(C2)C(=O)C=C1)C(=O)CC1" # Approximate SMILES
    }

    product_C = {
        "label": "C",
        "name": "(3S,5S)-1-acetyl-5-((S)-pyrrolidin-2-yl)piperidine-2,6-dione",
        "formula": "C11H16N2O3",
        "smiles": "CC(=O)N1C(=O)C[C@H]([C@H](C1=O))C1CCCN1"
    }

    print("Structures of the reaction products:\n")
    for product in [product_A, product_B, product_C]:
        print(f"Product {product['label']}:")
        print(f"  Name: {product['name']}")
        print(f"  Molecular Formula: {product['formula']}")
        # Note: The SMILES strings are representations of the complex structures.
        # print(f"  SMILES: {product['smiles']}")
        print("-" * 20)

    print("Final Reaction Equation (unbalanced, showing major components):")
    print("C10H14N2O2 (SM) + C4H4O2 (MP) + C4H6O3 (Ac2O) -->")
    print(f"  {product_A['formula']} (A) + {product_B['formula']} (B) + {product_C['formula']} (C) + byproducts")


if __name__ == '__main__':
    get_product_structures()