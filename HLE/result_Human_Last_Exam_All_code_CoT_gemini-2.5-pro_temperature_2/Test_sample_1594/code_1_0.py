def solve_pisa_schools_puzzle():
    """
    This function explains the shape formed by the two highschools
    in Pisa when seen from above.
    """
    school_1 = "Istituto Galilei-Pacinotti"
    school_2 = "Istituto Ulisse Dini"

    # Explain the observation from a satellite view
    print("This puzzle requires observing the layout of the two schools from a satellite map.")
    print(f"1. The building of the {school_1} forms a long, straight structure.")
    print(f"2. The building of the {school_2} stands right next to it, forming another long structure at a right angle.")
    
    # Formulate a symbolic equation to represent the final shape
    print("\nSymbolic 'equation' of the shapes:")
    building_1_shape = "part 1 (vertical wing)"
    building_2_shape = "part 2 (horizontal wing)"
    final_letter = "L"

    print(f"{building_1_shape} + {building_2_shape} = The letter '{final_letter}'")

solve_pisa_schools_puzzle()