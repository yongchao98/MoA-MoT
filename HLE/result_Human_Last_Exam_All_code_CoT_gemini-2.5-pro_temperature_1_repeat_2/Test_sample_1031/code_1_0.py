def solve_disease_spread():
    """
    Calculates the minimum number of initially infected sites to infect an n x n grid.
    """
    n = 14

    print(f"The grid size is n = {n}.")
    print("The spreading rule implies that any 2x2 block of healthy cells is stable.")
    print("Therefore, the initial set must contain at least one cell from every 2x2 subgrid.")
    print("-" * 30)

    # Lower Bound Calculation
    print("Step 1: Calculate the Lower Bound")
    # The number of 2x2 squares is (n-1)x(n-1)
    num_squares_dim = n - 1
    # We can find a set of disjoint 2x2 squares to establish a lower bound.
    # These are squares with top-left corners at (2i, 2j).
    # The number of possible indices for i (or j) is ceil((n-1)/2).
    num_disjoint_rows = (num_squares_dim + 1) // 2
    num_disjoint_cols = (num_squares_dim + 1) // 2
    
    lower_bound = num_disjoint_rows * num_disjoint_cols

    print(f"We can place disjoint 2x2 squares in a grid pattern.")
    print(f"Number of rows of disjoint squares = ceil(({n}-1)/2) = {num_disjoint_rows}")
    print(f"Number of columns of disjoint squares = ceil(({n}-1)/2) = {num_disjoint_cols}")
    print(f"The number of disjoint 2x2 squares is {num_disjoint_rows} * {num_disjoint_cols} = {lower_bound}.")
    print(f"This means we need at least {lower_bound} sites. So, Lower Bound = {lower_bound}.")
    print("-" * 30)

    # Upper Bound Calculation
    print("Step 2: Find a Valid Configuration (Upper Bound)")
    # We can infect all sites (i,j) where i and j are both odd.
    # For a grid of size n, there are n/2 odd indices.
    num_odd_rows = n // 2
    num_odd_cols = n // 2
    
    upper_bound = num_odd_rows * num_odd_cols
    
    print(f"Consider the set of sites (i,j) where i and j are both odd.")
    print(f"Number of odd row indices = {n} // 2 = {num_odd_rows}")
    print(f"Number of odd column indices = {n} // 2 = {num_odd_cols}")
    print(f"The total number of sites in this configuration is {num_odd_rows} * {num_odd_cols} = {upper_bound}.")
    print(f"This configuration hits every 2x2 square, so it is a valid solution. So, Upper Bound = {upper_bound}.")
    print("-" * 30)
    
    # Conclusion
    print("Step 3: Conclusion")
    final_answer = upper_bound
    print(f"The lower bound ({lower_bound}) and the upper bound ({upper_bound}) are equal.")
    print(f"Therefore, the minimum number of initially infected sites is {final_answer}.")

solve_disease_spread()
<<<49>>>