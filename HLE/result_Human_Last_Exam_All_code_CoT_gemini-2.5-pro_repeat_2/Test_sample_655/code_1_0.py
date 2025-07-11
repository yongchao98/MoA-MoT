import math

def solve_pioneer_probe_packing():
    """
    This function calculates the optimal packing of energy balls in a container
    with a surface area of at most 1050 cm^2, maximizing the total energy.
    """

    # --- Step 1: Determine the Optimal Container ---
    # Based on analysis, a cubic box provides the most efficient and practical packing.
    # For a cube with side L, Surface Area = 6 * L^2.
    # 6 * L^2 <= 1050  => L^2 <= 175 => L <= 13.22 cm.
    # Since dimensions must be multiples of 0.5 cm, the optimal side length is 13.0 cm.
    container_type = "box"
    side_length = 13.0
    container_desc = f"{container_type} {side_length}x{side_length}x{side_length}"
    surface_area = 6 * side_length**2

    # --- Step 2: Pack Large Balls (r=2cm, d=4cm) ---
    # We use a simple cubic lattice. The centers are placed on a grid with spacing 4.0 cm.
    # To respect the boundary, a center 'c' must be in [r, L-r]. For large balls, [2.0, 11.0].
    # Possible coordinates along one axis: 2.0, 6.0, 10.0. This gives 3 positions.
    positions_per_axis = 3
    num_large_balls = positions_per_axis * positions_per_axis * positions_per_axis
    e_large = 10  # MJ

    # --- Step 3: Pack Small Balls (r=1cm, d=2cm) in Gaps ---
    # The largest gaps in a simple cubic lattice are at the body-centers of the cells.
    # A cell formed by 8 large balls is a 4x4x4 cube. Its center is a potential site.
    # E.g., for a cell with corners at (2,2,2), (6,6,6), etc., the site is at (4,4,4).
    # We check if a small ball fits by ensuring its distance to the nearest large ball
    # center is at least r_large (2.0) + r_small (1.0) = 3.0 cm.
    # Distance from site (4,4,4) to large ball center (2,2,2) is sqrt((4-2)^2 + (4-2)^2 + (4-2)^2).
    dist_to_large_center = math.sqrt(2.0**2 + 2.0**2 + 2.0**2) # This is ~3.46 cm.
    # Since 3.46 > 3.0, a small ball fits.
    # The 3x3x3 lattice of large balls creates a 2x2x2 lattice of these interstitial sites.
    gaps_per_axis = positions_per_axis - 1
    num_small_balls = gaps_per_axis * gaps_per_axis * gaps_per_axis
    e_small = 1  # MJ

    # --- Step 4: Calculate Total Energy ---
    total_energy = (num_large_balls * e_large) + (num_small_balls * e_small)

    print("--- Pioneer Probe Energy Packing Optimization ---")
    print(f"\nOptimal container choice: {container_type.capitalize()}")
    print(f"Dimensions: {side_length} x {side_length} x {side_length} cm")
    print(f"Surface Area: {surface_area:.2f} cm^2")
    
    print("\nPacking result:")
    print(f"Number of 2-cm balls: {num_large_balls}")
    print(f"Number of 1-cm balls: {num_small_balls}")
    
    print("\nTotal potential energy calculation:")
    print(f"{num_large_balls} (2-cm balls) * {e_large} MJ/ball + {num_small_balls} (1-cm balls) * {e_small} MJ/ball = {total_energy} MJ")
    
    # --- Final Answer Formatting ---
    # Format: [C]a;b where C=description, a=num_small, b=num_large
    final_answer_string = f"[{container_desc}]{num_small_balls};{num_large_balls}"
    print(f"\n<<<{final_answer_string}>>>")


solve_pioneer_probe_packing()