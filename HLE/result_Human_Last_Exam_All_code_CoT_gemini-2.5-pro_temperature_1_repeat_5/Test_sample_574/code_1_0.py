def solve_ship_name_puzzle():
    """
    This script solves the word puzzle to find a Culture series ship name.
    """

    # 1. Define the words derived from the clues
    word1 = "SOBER"
    word2 = "COUNSEL"

    # 2. Define the source words for the explanation
    source_word1 = "BOLEROS"
    source_word2 = "COUNSELLORS"

    # 3. Form the final ship name
    ship_name = f"{word1} {word2}"

    # 4. Explain the solution step-by-step
    print(f"The word for 'sleeveless garments' is '{source_word1}'.")
    print(f"Using its letters, we can form the word: '{word1}'.\n")

    print(f"The word for 'experienced individuals who guide others' is '{source_word2}'.")
    print(f"Using its letters, we can form the word: '{word2}'.\n")

    print(f"Combining these two words gives the Culture series ship name: '{ship_name}'.\n")

    # 5. Create the "equation" as requested, showing each letter
    equation_part1 = " + ".join(list(word1))
    equation_part2 = " + ".join(list(word2))

    print("The final equation is:")
    print(f"{equation_part1} + {equation_part2} = {ship_name}")

# Execute the function
solve_ship_name_puzzle()
<<<SOBER COUNSEL>>>