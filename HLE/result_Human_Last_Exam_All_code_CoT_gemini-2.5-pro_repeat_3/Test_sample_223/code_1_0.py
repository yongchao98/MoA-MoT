def solve_chair_puzzle():
    """
    Calculates the maximum number of occupied chairs based on the puzzle's rules.
    """
    # Total number of chairs in the row.
    total_chairs = 20

    # The number of occupants increases only when a person sits in a chair with
    # no neighbors. This leads to a pattern where occupied chairs are separated
    # by at least one empty chair: O _ O _ ...

    # To find the maximum number of occupants, we need to find how many chairs can
    # be filled in this pattern. This is equivalent to calculating the ceiling
    # of the total number of chairs divided by 2.

    numerator = total_chairs
    denominator = 2

    # In Python, ceiling division can be performed using -(-a // b)
    max_occupied = -(-numerator // denominator)

    print("To find the maximum number of occupied chairs, we must find the densest packing")
    print("of people where no two are immediate neighbors. This corresponds to the pattern 'O_O_O_...'.")
    print("\nFor a row of 20 chairs, the calculation is the ceiling of the total chairs divided by 2.")
    print(f"The equation is: ceil({numerator} / {denominator})")
    print(f"The result is: {max_occupied}")

solve_chair_puzzle()