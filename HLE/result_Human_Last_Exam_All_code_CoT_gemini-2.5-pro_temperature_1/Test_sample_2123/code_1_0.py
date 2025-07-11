import sys

def solve_riddle():
    """
    This function solves the riddle by identifying a phonetic pun.
    """
    # The riddle's solution lies in the pun between the English word "wide"
    # (describing Vienna's boulevards) and a poet's surname.
    poet_surname = "Wilde"

    print("The riddle is based on a phonetic pun.")
    print("The wide boulevards of Vienna are described by the English word 'wide'.")
    print("This word sounds identical to the surname of the famous poet, Oscar Wilde.")
    print(f"\nThe surname is: {poet_surname}")

    # The prompt requests that the final code output each number in a final equation.
    # We will create a symbolic equation by enumerating the letters of the surname.
    print("\n---")
    print("To fulfill the 'equation' requirement, here is a breakdown of the name:")
    
    equation_parts = []
    # Create a part of the equation for each letter, e.g., "W(1)", "I(2)", etc.
    for i, char in enumerate(poet_surname):
        # We use i + 1 because enumeration starts at 0.
        num = i + 1
        equation_parts.append(f"{char}({num})")
    
    # Join the parts with a '+' to form the equation string.
    equation_str = " + ".join(equation_parts)
    
    # The result of the "equation" is the total number of letters.
    result = len(poet_surname)
    
    # Print the final equation, which includes each number (1, 2, 3, 4, 5).
    print(f"{equation_str} = {result} letters")

# Execute the function to solve the riddle.
solve_riddle()