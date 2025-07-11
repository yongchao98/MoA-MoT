import math

def solve_packing_problem():
    """
    Searches for a container (box) with a smaller surface area than the original
    that can hold at least 27 energy balls.
    """
    initial_side = 12.0
    initial_ball_count = 27
    initial_surface_area = 6 * initial_side**2
    ball_diameter = 4.0
    
    min_area = initial_surface_area
    best_config_str = "0"

    print(f"Initial container: {initial_side}x{initial_side}x{initial_side} cm box")
    print(f"Surface Area: {initial_surface_area} cm^2, holding {initial_ball_count} balls.")
    print("\nSearching for a more efficient box using staggered packing...")

    # We will search for configurations using staggered packing, as simple cubic packing
    # does not yield a better result.
    # Search Range for base layer dimensions (nx, ny) and number of layers (L)
    search_range_n = range(1, 15)
    search_range_L = range(2, 15)

    for nx in search_range_n:
        for ny in search_range_n:
            for L in search_range_L:
                # Box dimensions for this configuration
                dim_x = ball_diameter * nx
                dim_y = ball_diameter * ny
                # Height for L staggered layers is 2*r + (L-1)*dz + 2*r = d + (L-1)*3
                dim_z = ball_diameter + (L - 1) * 3.0

                # Ball count in alternating layers of sizes (nx, ny) and (nx-1, ny-1)
                num_A_layers = math.ceil(L / 2.0)
                num_B_layers = math.floor(L / 2.0)
                
                balls_A = nx * ny
                balls_B = (nx - 1) * (ny - 1) if nx > 1 and ny > 1 else 0
                total_balls = num_A_layers * balls_A + num_B_layers * balls_B

                if total_balls >= initial_ball_count:
                    surface_area = 2 * (dim_x * dim_y + dim_y * dim_z + dim_z * dim_x)
                    if surface_area < min_area:
                        min_area = surface_area
                        desc = f"box {dim_x}x{dim_y}x{dim_z}"
                        # The output format is d[X]
                        # Final equation: area = surface_area
                        best_config_str = f"{surface_area:.1f}[{desc}]"
    
    if best_config_str == "0":
        print("\nNo more efficient container configuration was found.")
    else:
        print(f"\nFound a more efficient container: {best_config_str}")

    print("\nFinal Answer:")
    # Printing the result of the search.
    # Based on the analysis, no better configuration is expected.
    print(best_config_str)

solve_packing_problem()