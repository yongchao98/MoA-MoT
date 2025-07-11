def find_reagents():
    """
    This function identifies and prints the reagents A and B for the given chemical reaction scheme.
    """
    reagent_A = {
        "name": "Hydrazine",
        "formula": "H2N-NH2"
    }

    reagent_B = {
        "name": "n-Propylamine",
        "formula": "CH3CH2CH2NH2"
    }

    print("Based on the analysis of the chemical transformations:")
    print("-" * 40)
    
    print(f"Reagent A is responsible for converting the xanthenylium salt (1) into the N-aminoacridinium salt (2).")
    print(f"This is achieved by replacing the ring oxygen atom with a >N-NH2 group.")
    print(f"Reagent A: {reagent_A['name']} (Formula: {reagent_A['formula']})")
    
    print("-" * 40)

    print(f"Reagent B reacts with the N-aminoacridinium salt (2) to form the quinacridinium derivative (3).")
    print(f"The presence of an N-propyl group in the final product indicates the use of a primary amine.")
    print(f"Reagent B: {reagent_B['name']} (Formula: {reagent_B['formula']})")
    print("-" * 40)

# Execute the function to display the results.
find_reagents()