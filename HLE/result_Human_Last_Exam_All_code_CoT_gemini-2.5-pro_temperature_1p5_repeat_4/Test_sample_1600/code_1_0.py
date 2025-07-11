def solve_feynman_diagrams():
    """
    Calculates a(3) for the number of non-vanishing Feynman diagrams in QED.
    
    The number of diagrams a(n) of order 2n follows the recurrence relation:
    n * a(n) = (4*n - 6) * a(n-1) - 4*(n - 3) * a(n-2)
    with initial values a(1) = 1 and a(2) = 1.
    """
    n = 3
    a1 = 1  # Value of a(n-2) when n=3
    a2 = 1  # Value of a(n-1) when n=3

    print(f"To find a({n}), we use the recurrence relation with n={n}:")
    print("n * a(n) = (4*n - 6) * a(n-1) - 4*(n - 3) * a(n-2)")
    print("\nWith a(1)=1 and a(2)=1, we can substitute n=3:")
    
    # Building the equation string step-by-step
    c1 = 4 * n - 6
    c2 = 4 * (n - 3)
    
    print(f"{n} * a({n}) = (4*{n} - 6) * a(2) - 4*({n} - 3) * a(1)")
    print(f"{n} * a({n}) = ({c1}) * {a2} - ({c2}) * {a1}")
    
    term1_val = c1 * a2
    term2_val = c2 * a1
    numerator = term1_val - term2_val
    
    print(f"{n} * a({n}) = {term1_val} - {term2_val}")
    print(f"{n} * a({n}) = {numerator}")
    
    result = numerator / n
    print(f"a({n}) = {numerator} / {n}")
    print(f"a({n}) = {int(result)}")

solve_feynman_diagrams()