def identify_reagents():
    """
    This function analyzes the provided chemical reaction scheme and identifies
    the reagents A and B based on the structural transformations shown.
    """
    # Analysis of the transformation from compound 1 to 2
    reagent_A_name = "Hydrazine"
    reagent_A_formula = "N2H4"
    
    # Analysis of the transformation from compound 2 to 3
    reagent_B_name = "n-Propylamine"
    reagent_B_formula = "CH3CH2CH2NH2"
    
    # Print the results
    print("Based on the chemical transformations shown in the diagram:")
    print("-" * 50)
    
    # Details for Reagent A
    print("To convert compound 1 to compound 2:")
    print("The reaction involves replacing a C-O-C bridge with a C-N(NH2)-C bridge.")
    print(f"Reagent A is: {reagent_A_name} ({reagent_A_formula})")
    print("-" * 50)

    # Details for Reagent B
    print("To convert compound 2 to compound 3:")
    print("The reaction involves replacing a C-O-C bridge with a C-N(propyl)-C bridge and a deamination of the N-NH2 group.")
    print(f"Reagent B is: {reagent_B_name} ({reagent_B_formula})")
    print("-" * 50)

# Execute the function to display the answer
identify_reagents()