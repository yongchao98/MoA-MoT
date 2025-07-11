def find_dinosaur_hydrocarbon():
    """
    This function identifies and prints the name of a polycyclic compound
    named after a prehistoric creature.
    """
    
    # The creature that went extinct over 65 million years ago
    creature_name = "Pterodactyl"
    
    # The name of the polycyclic compound
    compound_name = "Pterodactyladiene"
    
    # The chemical formula of the compound
    # It contains Carbon (C), Hydrogen (H), and Nitrogen (N) atoms.
    chemical_formula = {
        "element_C": 13,
        "element_H": 21,
        "element_N": 1
    }

    print(f"The polycyclic compound named after a creature that has been extinct for over 65 million years is {compound_name}.")
    print(f"It is named after the {creature_name}, as its chemical structure is said to resemble one.")
    print("\nIts chemical formula is C13H21N, which breaks down as follows:")
    
    # Printing the numbers from the 'equation' (chemical formula)
    print(f"Carbon atoms: {chemical_formula['element_C']}")
    print(f"Hydrogen atoms: {chemical_formula['element_H']}")
    print(f"Nitrogen atoms: {chemical_formula['element_N']}")

# Execute the function to display the answer
find_dinosaur_hydrocarbon()