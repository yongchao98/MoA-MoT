def solve_vienna_riddle():
    """
    This function solves the riddle by combining descriptive elements of Vienna
    to form the surname of an English poet.
    """
    # The riddle points to the characteristic color of the buildings on Vienna's main boulevard.
    descriptive_color = "Brown"

    # The main boulevard in Vienna is the Ringstrasse, or "Ring" for short.
    # The "-ing" part of the poet's name comes from "Ring".
    boulevard_name_part = "ing"

    # The surname is formed by combining the color and the part from the boulevard's name.
    poet_surname = descriptive_color + boulevard_name_part

    print("The riddle describes the boulevards of Vienna.")
    print("The logic combines two key observations:")
    print(f"1. The color of the imperial buildings: '{descriptive_color}'")
    print(f"2. The name of the main boulevard, the 'Ring', which gives us: '{boulevard_name_part}'")
    print("-" * 20)
    print(f"Combining these gives the surname: {descriptive_color} + {boulevard_name_part} = {poet_surname}")
    print("-" * 20)
    print(f"The English poet is Robert {poet_surname}.")

solve_vienna_riddle()