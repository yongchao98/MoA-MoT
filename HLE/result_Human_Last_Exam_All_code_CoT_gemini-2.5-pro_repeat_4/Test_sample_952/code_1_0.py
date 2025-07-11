import math

def solve_diamond_problem():
    """
    This function calculates the largest value k based on the problem description.
    """
    # The side length of the square table.
    n = 2024

    # The problem of finding a "stuck" arrangement of diamonds is equivalent to
    # finding a "maximal independent set" (MIS) in the grid graph.
    # We are looking for the largest k such that any arrangement of k diamonds
    # is NOT a stuck arrangement (i.e., not an MIS).
    # This k is one less than the size of the smallest possible MIS.

    # For an n x n grid graph (with n >= 9), the minimum size of a maximal
    # independent set is given by the formula: ceil((n*n) / 5).
    
    # 1. Calculate the total number of cells in the grid.
    total_cells = n * n
    
    # 2. Calculate the minimum size of a maximal independent set (M_min).
    # This is the size of the smallest "stuck" arrangement.
    M_min = math.ceil(total_cells / 5)
    
    # 3. The largest value k for which any arrangement is guaranteed to have a
    # mobile diamond is M_min - 1.
    k = M_min - 1
    
    # Output the steps of the calculation as requested.
    # The final equation is k = ceil(n*n / 5) - 1.
    # Here we show the numbers involved in this equation.
    print(f"Grid size: {n} x {n}")
    print(f"Total cells: {total_cells}")
    print(f"Formula for the minimum size of a 'stuck' arrangement (M_min): ceil(Total Cells / 5)")
    print(f"M_min = ceil({total_cells} / 5) = {M_min}")
    print(f"The value k is M_min - 1.")
    print(f"Final equation: {M_min} - 1 = {k}")
    
    # The final answer
    print(f"\nThe largest value is {k}.")


solve_diamond_problem()