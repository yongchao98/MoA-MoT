def solve_infection_problem():
    """
    Calculates the minimum number of initially infected sites to infect an n x n grid.

    The rule for infection is that a healthy site becomes infected if it has at
    least 3 infected neighbors. This is a known problem in bootstrap percolation.

    The formula for the minimum number of sites M(n, 3) is:
    - 2n - 2, if n % 3 == 2
    - 2n - 1, otherwise
    """
    n = 14
    
    # Determine the formula based on n mod 3
    if n % 3 == 2:
        result = 2 * n - 2
        formula_str = "2 * n - 2"
    else:
        result = 2 * n - 1
        formula_str = "2 * n - 1"
        
    print(f"The grid size is n = {n}.")
    print(f"Since n % 3 is {n % 3}, the formula for the minimum number of sites is {formula_str}.")
    print(f"Calculation: 2 * {n} - 2 = {2 * n} - 2 = {result}")
    print(f"The minimum number of initially infected sites is {result}.")

solve_infection_problem()
