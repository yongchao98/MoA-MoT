def identify_compound():
    """
    This script identifies the provided alkaloid and displays its properties.
    """
    
    # The compound is identified as Matrine.
    # Its chemical formula is C15H24N2O.
    
    name = "Matrine"
    atom_counts = {
        "Carbon": 15,
        "Hydrogen": 24,
        "Nitrogen": 2,
        "Oxygen": 1,
    }

    print("Analyzing the molecular formula: C15H24N2O")
    print("The count of each atom is as follows:")
    
    # Printing each number from the "equation" (molecular formula)
    for atom, count in atom_counts.items():
        print(f"Number of {atom} atoms: {count}")
    
    print("\n--------------------------------")
    print("The name of the alkaloid is:")
    print(name)
    print("--------------------------------")

identify_compound()