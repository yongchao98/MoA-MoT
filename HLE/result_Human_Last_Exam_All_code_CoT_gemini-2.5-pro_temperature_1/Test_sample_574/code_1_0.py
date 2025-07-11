def solve_culture_ship_puzzle():
    """
    This script solves the word puzzle to reveal a Culture series ship name.
    """

    # Clue 1: "experienced and trusted individuals who guide and advise others"
    source_word_1 = "MENTORS"
    derived_word_1 = "STERN"

    # Clue 2: "sleeveless garments that drape over the back and shoulders"
    source_word_2 = "ACADEMIC ROBES"
    derived_word_2 = "CANDOR"

    # The final ship name is the combination of the two derived words.
    ship_name_part_1 = derived_word_1
    ship_name_part_2 = derived_word_2

    print("Solving the Culture ship puzzle:")
    print("-" * 35)
    print(f"The word for 'experienced individuals' is '{source_word_1}'.")
    print(f"From its letters, we can make the word: '{derived_word_1}'")
    print("-" * 35)
    print(f"The word for 'sleeveless garments' is '{source_word_2}'.")
    print(f"From its letters, we can make the word: '{derived_word_2}'")
    print("-" * 35)
    print("Combining these two words gives the ship name.")
    # The prompt asked to output each part of the final "equation" or name.
    print(f"The final ship name is: {ship_name_part_1} {ship_name_part_2}")

solve_culture_ship_puzzle()