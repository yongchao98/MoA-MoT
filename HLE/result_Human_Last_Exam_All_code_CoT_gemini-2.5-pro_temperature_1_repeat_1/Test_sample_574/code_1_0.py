def solve_riddle():
    """
    This function solves the word riddle and prints the steps and the final answer.
    """
    # Step 1: Identify the source words from the clues.
    garment_word = "CAPES"
    advisor_word = "MENTORS"

    # Step 2: Form new words using letters from the source words.
    first_word = "SPACE"
    second_word = "STORM"

    # Step 3: Combine the new words to get the final ship name.
    ship_name = f"{first_word} {second_word}"

    # Step 4: Print the logic and the final answer.
    print(f"The word for 'sleeveless garments' is '{garment_word}'. From its letters, we can make the word '{first_word}'.")
    print(f"The word for 'trusted individuals' is '{advisor_word}'. From its letters, we can make the word '{second_word}'.")
    print("\nCombining these gives a ship name in the style of the Culture series.")
    
    # Per the instructions, output the final equation showing each component.
    print(f"\nFinal Equation:")
    print(f"{first_word} + {second_word} = {ship_name}")

solve_riddle()