import math

def solve_cutting_problem():
    """
    Solves the billet cutting problem based on the user's formulation.
    """
    
    # --- Problem Definition ---
    
    # Billet dimensions in cm and grid units (1 grid unit = 0.5 cm)
    length_cm, width_cm, height_cm = 16, 11, 4
    length_grid = length_cm * 2  # 32
    width_grid = width_cm * 2    # 22
    height_grid = height_cm * 2  # 8

    # T1 Cube properties
    price_t1 = 5
    # Valid range for T1 centers
    x_min_t1, x_max_t1 = 1, 31
    y_min_t1, y_max_t1 = 1, 21
    z_min_t1, z_max_t1 = 1, 7
    # Separation distance from the user's T1-T1 constraint: min(dx,dy,dz) >= 2
    # This implies dx >= 2, dy >= 2, and dz >= 2.
    separation_t1 = 2

    # --- Analysis ---
    
    print("Based on the provided (and non-standard) problem formulation, the optimal strategy is to fill the billet with T1 cubes.")
    print("The constraint 'min(|xi-xj|, |yi-yj|, |zi-zj|) >= 2' for T1-T1 allows for a simple grid-based packing.")
    print("-" * 20)

    # --- Calculation for T1-only strategy ---

    # Calculate the number of possible center positions along each axis
    # Formula: floor((max_coord - min_coord) / separation) + 1
    num_x_t1 = (x_max_t1 - x_min_t1) // separation_t1 + 1
    num_y_t1 = (y_max_t1 - y_min_t1) // separation_t1 + 1
    num_z_t1 = (z_max_t1 - z_min_t1) // separation_t1 + 1
    
    # Total number of T1 cubes is the product of positions along each axis
    total_t1 = num_x_t1 * num_y_t1 * num_z_t1
    
    # Total value from this strategy
    max_value = total_t1 * price_t1

    # --- Output Results ---
    
    print(f"Number of T1 cubes that can be placed along the X-axis: {num_x_t1}")
    print(f"Number of T1 cubes that can be placed along the Y-axis: {num_y_t1}")
    print(f"Number of T1 cubes that can be placed along the Z-axis: {num_z_t1}")
    print("-" * 20)
    
    # The final equation as requested by the user
    print("The highest value is achieved with the following calculation:")
    print(f"Number of T1 Cubes * Price per T1 Cube = Total Value")
    print(f"{total_t1} * {price_t1} = {max_value}")

solve_cutting_problem()
