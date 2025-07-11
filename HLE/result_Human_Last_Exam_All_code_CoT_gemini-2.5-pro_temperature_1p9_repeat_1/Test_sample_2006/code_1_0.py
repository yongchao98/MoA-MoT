def solve_riddle():
    """
    Solves the base equation derived from the haiku riddle.

    The riddle "8 + 30 = 41" is translated into the algebraic equation:
    8 + (3*b + 0) = (4*b + 1)
    This simplifies to: 8 + 3b = 4b + 1
    We can solve for b programmatically.
    """

    # We can solve this algebraically (b = 7), but let's use a loop
    # to demonstrate a coding-based solution. The base must be larger
    # than any single digit used in the numbers (i.e., b > 4).
    found_base = None
    for b in range(5, 100):
        # Represent the value of "30" in base b
        val_30 = 3 * b + 0
        # Represent the value of "41" in base b
        val_41 = 4 * b + 1
        
        # Check if the equation 8 + val_30 == val_41 holds true
        if 8 + val_30 == val_41:
            found_base = b
            break

    if found_base is not None:
        num1 = 8
        num2 = 30
        result = 41
        print("The riddle describes an equation in an unknown base.")
        print(f"The equation is: {num1} + {num2} = {result}")
        print(f"Solving for the base reveals that the equation is true in base {found_base}.")
        print(f"Final equation: {num1} + {num2} = {result} (in base {found_base})")
    else:
        print("Could not find a base that solves the riddle.")

solve_riddle()