def solve_riddle():
    """
    This function solves the riddle by breaking it down and explaining the logic.
    """
    
    # The two words that form the solution
    word_x = "Felt"
    word_y = "Hat"

    print("Solving the riddle by breaking it into parts:")
    print("-" * 40)
    
    # Explanation for the first word
    print("Clue 1: An event on June 20, 2019, shows Nikita Kucherov winning a 'hat trick' of 3 NHL awards.")
    print(f"From the emotions of the winner, we can derive the first word X: '{word_x}'")
    print("\n")

    # Explanation for the second word
    print("Clue 2: An object whose manufacturing used toxic mercury salts.")
    print("This refers to old-fashioned hats, whose felt was treated with mercury.")
    print(f"From the 'hat trick' and the historical clue, we derive the second word Y: '{word_y}'")
    print("-" * 40)
    
    # Final answer composition
    print("Combining the two words gives the final answer.")
    print(f"The final answer is: {word_x} {word_y}")

solve_riddle()