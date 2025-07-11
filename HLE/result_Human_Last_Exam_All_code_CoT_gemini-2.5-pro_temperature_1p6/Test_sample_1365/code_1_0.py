def solve():
    """
    This function calculates and prints the upper bound on the number of mistakes
    made by the majority voting algorithm based on the derived formula.
    """
    
    # The variables n (number of experts) and c (mistake threshold) are symbolic.
    # The derived formula provides an upper bound M on the number of algorithm mistakes.
    # M = M1 + M2
    # M1 <= c - 1
    # M2 <= (n - 1) * c / 2
    # So, M <= (c - 1) + (n - 1) * c / 2

    c = 'c'
    n = 'n'
    one = 1
    two = 2
    
    # Print the derived formula for the upper bound.
    # The format displays each component of the final equation.
    print("An upper bound M on the number of mistakes is given by:")
    print(f"M <= ({c} - {one}) + ({c} * ({n} - {one})) / {two}")

solve()