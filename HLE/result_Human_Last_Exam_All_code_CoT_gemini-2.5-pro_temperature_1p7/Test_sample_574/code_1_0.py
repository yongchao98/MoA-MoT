def solve_culture_ship_puzzle():
    """
    This function solves the word puzzle and prints the step-by-step solution.
    """

    # Step 1: Identify the source phrase and the first word.
    word1 = "Sleeper"
    source_phrase1 = "REVERSIBLE CAPE"
    description1 = "sleeveless garments that drape over the back and shoulders"

    print(f"Step 1: The phrase describing '{description1}' is '{source_phrase1}'.")
    print(f"Using letters from '{source_phrase1}', we can form the word: '{word1}'.")
    print("-" * 20)

    # Step 2: Identify the source phrase and the second word.
    word2 = "Service"
    source_phrase2 = "SPECIAL ADVISER"
    description2 = "experienced and trusted individuals who guide and advise others"

    print(f"Step 2: The phrase for '{description2}' is '{source_phrase2}'.")
    print(f"Using letters from '{source_phrase2}', we can form the word: '{word2}'.")
    print("-" * 20)

    # Step 3: Combine the words to get the final ship name.
    final_ship_name = f"{word1} {word2}"
    print("Step 3: Combine the two words to form the final ship name.")
    
    # Per the instructions, printing each component of the final equation.
    print(f"The final equation is: {word1} + {word2} = {final_ship_name}")

solve_culture_ship_puzzle()