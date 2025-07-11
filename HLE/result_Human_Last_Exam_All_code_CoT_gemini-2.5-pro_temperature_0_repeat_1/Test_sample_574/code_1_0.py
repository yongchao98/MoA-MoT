def solve_culture_ship_riddle():
    """
    This function solves the word puzzle and prints the result.
    """

    # The first word is formed from the letters of "PEERLESS",
    # which describes "experienced and trusted individuals who guide and advise others".
    word_from_advisors = "SLEEPER"

    # The second word is formed from the letters of "VICARESSES",
    # describing the wearers of "sleeveless garments that drape over the back and shoulders".
    word_from_garments = "SERVICE"

    # The two words are combined to form the ship name, starting with 'S'.
    # The instruction "output each number in the final equation" is interpreted
    # as showing each component word in the final result.
    final_ship_name = f"{word_from_advisors} {word_from_garments}"

    print(f"{word_from_advisors} + {word_from_garments} = {final_ship_name}")

solve_culture_ship_riddle()