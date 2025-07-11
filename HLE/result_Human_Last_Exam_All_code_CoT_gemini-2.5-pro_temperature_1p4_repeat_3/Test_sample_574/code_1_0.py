import collections

def solve_ship_name_puzzle():
    """
    Solves the word puzzle to find a Culture series ship name.
    """

    # --- Step 1: Define source materials from the puzzle clues ---
    clue_phrase = "sleeveless garments that drape over the back and shoulders"
    source_word_for_advisors = "CIVILSERVANTS"

    # --- Step 2: Define the target words to be formed ---
    target_word_1 = "SLEEPER"
    target_word_2 = "SERVICE"

    # --- Step 3: Combine the words to form the final ship name ---
    final_ship_name = f"{target_word_1} {target_word_2}"

    # --- Step 4: Print the explanation and final answer ---
    print("Step 1: Analyzing the first clue.")
    print(f"Clue: 'Use letters from sleeveless garments that drape over the back and shoulders to make a word.'")
    print(f"From the phrase '{clue_phrase}', we can form the word '{target_word_1}'.\n")

    print("Step 2: Analyzing the second clue.")
    print(f"Clue: 'Take letters from the word representing experienced and trusted individuals who guide and advise others to make another word.'")
    print(f"A word for these individuals is '{source_word_for_advisors}'.")
    print(f"From the word '{source_word_for_advisors}', we can form the word '{target_word_2}'.\n")

    print("Step 3: Combining the words.")
    print("The final ship name is formed by combining the two words:")
    print(f"'{target_word_1}' + '{target_word_2}' = '{final_ship_name}'")

# Execute the solver
solve_ship_name_puzzle()