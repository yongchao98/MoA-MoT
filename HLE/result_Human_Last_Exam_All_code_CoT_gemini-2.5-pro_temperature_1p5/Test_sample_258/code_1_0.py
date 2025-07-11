import math

def solve_grid_circle_problem():
    """
    Calculates the minimal and maximal numbers of grid cells a circle
    of radius 500 can cross.
    """
    R = 500

    # 1. Calculate the maximal number of cells
    # This occurs for a 'generic' placement, C_max = 8R
    max_cells = 8 * R
    
    # 2. Calculate the minimal number of cells
    # We need to find k, the number of unique pairs of positive integers {m, n}
    # such that m^2 + n^2 = R^2. C_min = 8R - 4k
    R2 = R * R
    pairs = set()
    
    # We only need to check for m from 1 up to floor(R/sqrt(2))
    # to find unique pairs {m,n} where m <= n.
    limit = int(R / math.sqrt(2))
    
    for m in range(1, limit + 1):
        n2 = R2 - m * m
        # Check if n2 is a perfect square
        if n2 > 0:
            n = int(math.sqrt(n2))
            if n * n == n2:
                # We found a pair of positive integers (m, n)
                # Since we iterate m up to R/sqrt(2), m will be <= n.
                # Adding the sorted tuple to a set handles uniqueness.
                pairs.add(tuple(sorted((m, n))))
    
    k = len(pairs)
    min_cells = 8 * R - 4 * k
    
    # 3. Print the results clearly, including the components of the equations.
    print(f"The radius of the circle is R = {R}.")
    
    print("\n--- Maximal Number of Cells ---")
    print(f"The maximal number of cells is given by the formula C_max = 8 * R.")
    print(f"C_max = 8 * {R} = {max_cells}")
    
    print("\n--- Minimal Number of Cells ---")
    print("The minimal number of cells is given by C_min = 8 * R - 4 * k, where k is the number of unique pairs of positive integers {m, n} such that m^2 + n^2 = R^2.")
    print(f"For R = {R}, we need to solve m^2 + n^2 = {R}^2 = {R2}.")
    
    # Python code found k pairs. Let's list them.
    sorted_pairs = sorted(list(pairs))
    print(f"The {k} pairs are: {sorted_pairs}")
    
    print(f"So, k = {k}.")
    print(f"C_min = 8 * {R} - 4 * {k} = {8 * R} - {4 * k} = {min_cells}")
    
    print("\n--- Final Answer ---")
    print(f"The minimal and maximal numbers of grid's cells are ({min_cells}, {max_cells}).")

solve_grid_circle_problem()
<<< (3988, 4000) >>>