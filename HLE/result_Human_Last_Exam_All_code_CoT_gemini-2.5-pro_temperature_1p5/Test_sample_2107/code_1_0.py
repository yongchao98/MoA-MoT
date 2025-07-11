def solve_riddle():
    """
    This function solves the multi-step riddle to find the value of X12.
    """
    # Step 1: Define the values for X1 through X9 based on the clues.
    X1 = 5
    X2 = 2
    X3 = 7
    X4 = 8
    X5 = 10
    X6 = 10
    X7 = 2
    X8 = 1  # This value makes the first equation balance.
    X9 = 4

    # Step 2: Verify the first equation.
    # X1's X2 + (X3 X4 - X5's X6) + X7's X4 = X8's X9
    # 's is multiplication, and ' ' is string concatenation.
    lhs = X1 * X2 + (int(str(X3) + str(X4)) - X5 * X6) + X7 * X4
    rhs = X8 * X9
    
    # This verification should result in True.
    # print(f"Verification that the first equation holds: {lhs} = {rhs} is {lhs == rhs}")

    # Step 3: Define the values for X10 and X11 based on word lengths.
    # X10: From Prutkov's quote ending in "посторонний" (stranger).
    X10_word = "посторонний"
    X10 = len(X10_word)

    # X11: From Belza's "bibliography".
    X11_word = "библиография"
    X11 = len(X11_word)
    
    # Step 4: Calculate X12 from the second phrase.
    # X8's X9's X10 X11 X12.
    # We interpret this as X12 being the sum of the previous terms.
    X12 = X8 + X9 + X10 + X11
    
    # Final output as requested
    print("The final equation to find X12 is a sum of the derived numbers:")
    print(f"{X8} + {X9} + {X10} + {X11} = {X12}")
    
    # The final answer for X12 is the result of the sum.
    # No extra text, just the return value as a string per instructions.
    # But the puzzle also asks to return the answer using a special format.
    # The format is <<<answer>>>, so we will return it that way.

solve_riddle()
print("<<<28>>>")