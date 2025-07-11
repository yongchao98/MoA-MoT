#
# This script solves a Chinese character riddle by analyzing its clues.
#

def solve_riddle():
    """
    Analyzes the riddle's clues and prints the solution.
    """
    print("Solving the Chinese character riddle step-by-step:")
    print("-" * 40)

    # Explain the riddle's logic
    print("The riddle is a classic example of wordplay.")
    print("While it seems to describe a character with 3 horizontal and 3 vertical strokes, this is a trick.")
    print("\nLet's break it down:")
    
    # Clue 1: Horizontal Strokes
    print("\n1. 'One horizontal stroke, another horizontal stroke...': This describes the two horizontal lines (二) in the character.")
    
    # Clue 2: Vertical Strokes
    print("2. 'one vertical on the left, one vertical on the right...': This perfectly describes the two vertical lines (||) that frame the character.")

    # The Trick
    print("\n3. The '...after another' part is the trick! The riddle-maker adds an extra stroke in the description to mislead you.")

    # The 'Equation'
    print("\nThe visual 'equation' for the character is:")
    print("Number of horizontal strokes = 2")
    print("Number of vertical strokes   = 2")
    print("\nThese strokes are arranged to form a grid shape.")

    # The final character
    final_char = "井"
    print("\nWhen you combine two horizontal strokes with two vertical strokes ('one on the left, one on the right'), you get the character:")
    print(f"'{final_char}' (jǐng), which means 'well'.")

# Execute the function to print the solution
solve_riddle()