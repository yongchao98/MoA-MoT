def identify_alkaloid():
    """
    This function identifies the alkaloid compound from the provided structure
    and prints its name and chemical formula.
    """
    compound_name = "Matrine"
    # The chemical formula is C15H24N2O
    atoms = {
        "Carbon": 15,
        "Hydrogen": 24,
        "Nitrogen": 2,
        "Oxygen": 1
    }

    print(f"The name of the alkaloid compound is: {compound_name}")
    print("\nThe chemical formula is C15H24N2O.")
    print("The count for each element is as follows:")
    
    # Printing each number from the formula as requested
    for element, count in atoms.items():
        print(f"- {element}: {count}")

identify_alkaloid()