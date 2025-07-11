def solve_character_riddle():
    """
    This function analyzes the riddle and prints a step-by-step explanation
    to identify the Chinese character.
    """
    print("Let's solve the riddle by analyzing its parts:")
    print("-" * 30)

    # Clue 1: Horizontal strokes
    clue1 = "One horizontal stroke, another horizontal stroke, after another"
    analysis1 = "This part describes multiple horizontal strokes, similar to those in the character '三' (sān, meaning three)."
    print(f"Clue 1: '{clue1}'")
    print(f"Analysis: {analysis1}\n")

    # Clue 2: Vertical strokes and their positions
    clue2 = "one vertical on the left, one vertical on the right"
    analysis2 = "This is a key structural clue. It clearly indicates that the character has two main vertical strokes forming its left and right boundaries."
    print(f"Clue 2: '{clue2}'")
    print(f"Analysis: {analysis2}\n")

    # Clue 3: Repetition
    clue3 = "one vertical stroke, another vertical stroke, after another"
    analysis3 = "The repetition in the riddle emphasizes the presence of multiple strokes and the character's symmetry."
    print(f"Clue 3 & Repetition: '{clue3}'")
    print(f"Analysis: {analysis3}\n")
    
    # Conclusion
    print("-" * 30)
    print("Conclusion:")
    print("Combining the clues, we are looking for a character with a vertical stroke on the left and another on the right, with multiple horizontal strokes connected to them.")
    print("The character that perfectly matches this description is '非' (fēi).")
    print("It is composed of a left part (a vertical line with three horizontal lines) and a right part (also a vertical line with three horizontal lines).")
    print("\nThe character described by the riddle is:")
    print("非")

# Execute the function to solve the riddle.
solve_character_riddle()