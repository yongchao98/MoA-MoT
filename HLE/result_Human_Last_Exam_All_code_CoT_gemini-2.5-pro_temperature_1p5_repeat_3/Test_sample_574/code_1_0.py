def solve_ship_name_puzzle():
    """
    Solves a word puzzle to find a Culture series ship name.
    """

    # The first word is derived from the letters of "CAPES".
    word_from_capes = "SCAPE"

    # The second word is derived from the letters of "MENTORS".
    word_from_mentors = "MENTOR"

    # The final ship name is the combination of the two words.
    final_ship_name = word_from_capes + word_from_mentors

    # Print the solution step-by-step.
    print("The first word is formed from the letters of 'CAPES' (sleeveless garments).")
    print(f"Word 1: {word_from_capes}")
    print("\nThe second word is formed from the letters of 'MENTORS' (trusted individuals).")
    print(f"Word 2: {word_from_mentors}")
    print("\nCombining these forms the Culture ship name.")
    print("The final equation is:")
    
    # Output each component of the final equation.
    print(f"{word_from_capes} + {word_from_mentors} = {final_ship_name}")

# Execute the function to print the result.
solve_ship_name_puzzle()