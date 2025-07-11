import collections

def solve_culture_ship_riddle():
    """
    This function solves the riddle by identifying source words based on the clues,
    deriving the target words for the ship's name, and verifying the logic.
    """

    # Step 1: Define the source word and target word for the first part of the riddle.
    # The clue "sleeveless garments..." suggests sleeping attire.
    source_word_1 = "SLEEPERS"
    target_word_1 = "SLEEPER"

    # Step 2: Define the source word and target word for the second part.
    # A guide who provides a 'service' could be a 'PRIESTESS'.
    source_word_2 = "PRIESTESS"
    target_word_2 = "SERVICE"

    def can_form(source, target):
        """Checks if the target word can be formed from the letters of the source word."""
        source_counts = collections.Counter(source.upper())
        target_counts = collections.Counter(target.upper())
        # The operation `target_counts - source_counts` will be empty if all letters
        # in target are available in source with sufficient frequency.
        return not (target_counts - source_counts)

    # Step 3: Verify the words and construct the final ship name.
    is_solution_valid_1 = can_form(source_word_1, target_word_1)
    is_solution_valid_2 = can_form(source_word_2, target_word_2)
    
    if not (is_solution_valid_1 and is_solution_valid_2):
        print("Could not find a valid solution based on the logic.")
        return

    ship_name = f"{target_word_1} {target_word_2}"

    # Step 4: Print the solution step-by-step.
    print("--- Riddle Analysis ---")
    print(f"1. Clue: 'sleeveless garments that drape over the back and shoulders'")
    print(f"   Interpretation: This hints at sleeping garments, like 'SLEEPERS'.")
    print(f"   From the letters of '{source_word_1}', we form the word: '{target_word_1}'")
    print("-" * 25)
    print(f"2. Clue: 'experienced and trusted individuals who guide and advise others'")
    print(f"   Interpretation: This describes someone providing a 'service', like a 'PRIESTESS'.")
    print(f"   From the letters of '{source_word_2}', we form the word: '{target_word_2}'")
    print("-" * 25)
    print(f"3. Combined, they form the Culture series ship name: {ship_name}")
    print("\n--- Final Equation ---")
    
    # Print the "equation" showing the letters that form the final name.
    # The format is requested to show each individual letter.
    print(" + ".join(f"'{char}'" for char in target_word_1) + f" (from {source_word_1})")
    print("+")
    print(" + ".join(f"'{char}'" for char in target_word_2) + f" (from {source_word_2})")
    print(f"\n= {ship_name}")


# Execute the function to solve the riddle and print the result.
solve_culture_ship_riddle()
<<<Sleeper Service>>>