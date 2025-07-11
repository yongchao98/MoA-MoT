def solve_ship_name_puzzle():
    """
    This function solves the word puzzle and prints the solution.
    """
    # The first clue refers to "ROBES".
    # The word derived from its letters is "Sober".
    first_word = "Sober"

    # The second clue refers to "COUNSELORS".
    # The word derived from its letters is "Counsel".
    second_word = "Counsel"

    # The final ship name is the combination of the two words.
    ship_name = f"{first_word} {second_word}"

    print("The first word is formed from the letters of 'ROBES' (garments):")
    print(first_word)
    print("\nThe second word is formed from the letters of 'COUNSELORS' (advisors):")
    print(second_word)
    print("\nThe final Culture ship name is:")
    # The prompt asks to output each part of the final equation.
    print(f"{first_word} + {second_word} = {ship_name}")

solve_ship_name_puzzle()