import math

def solve_pioneer_probe_packing():
    """
    Analyzes the Pioneer probe packing problem to find a more material-efficient container.
    """
    # --- Step 1: Analyze the Initial Box ---
    ball_radius_cm = 2.0
    ball_diameter_cm = 2.0 * ball_radius_cm
    initial_dim_cm = 12.0

    # For a dimension of 12 cm, we can place ball centers at 2, 6, and 10 cm.
    # This allows for 3 balls along each axis.
    # Capacity = 3 * 3 * 3 = 27 balls.
    initial_ball_capacity = 27
    initial_surface_area_cm2 = 6 * (initial_dim_cm ** 2)

    print("--- Initial Box Analysis ---")
    print(f"Original Box Dimensions: {int(initial_dim_cm)}x{int(initial_dim_cm)}x{int(initial_dim_cm)} cm")
    print(f"Surface Area: 2 * (12*12 + 12*12 + 12*12) = {int(initial_surface_area_cm2)} cm^2")
    print(f"Ball Capacity (Grid Packing): 3 x 3 x 3 = {initial_ball_capacity} balls")
    print("-" * 30)
    print("\n--- Searching for a More Efficient Box ---")

    # --- Step 2 & 3: Optimization Search ---
    best_solution = None
    min_found_surface_area = initial_surface_area_cm2

    # The search range for the number of balls per dimension.
    # A very elongated box is inefficient, so we can limit the search reasonably.
    search_limit = initial_ball_capacity + 1  # Search up to 27 balls in one line

    # Iterate through possible numbers of balls along each dimension (n_l, n_w, n_h)
    for n_l in range(1, search_limit):
        for n_w in range(n_l, search_limit):  # Start from n_l to avoid duplicate shapes (e.g., 8x4x2 vs 4x8x2)
            for n_h in range(n_w, search_limit):
                total_balls = n_l * n_w * n_h

                # Constraint 1: Must hold at least as many balls as the original box
                if total_balls >= initial_ball_capacity:
                    # Minimum dimension required for n balls in a line is 4*n
                    L = 4 * n_l
                    W = 4 * n_w
                    H = 4 * n_h

                    current_surface_area = 2 * (L*W + L*H + W*H)

                    # Constraint 2: Surface area must be smaller
                    if current_surface_area < min_found_surface_area:
                        min_found_surface_area = current_surface_area
                        best_solution = {
                            "dims": (L, W, H),
                            "area": current_surface_area
                        }

    # --- Step 4: Evaluate and Conclude ---
    if best_solution:
        # A solution was found
        dims = best_solution["dims"]
        area = int(best_solution["area"])
        a, b, c = sorted(dims)
        
        print(f"\nFound a more efficient box design!")
        print(f"Dimensions: {a}cm x {b}cm x {c}cm")
        print(f"Surface Area: 2 * ({a}*{b} + {a}*{c} + {b}*{c}) = {area} cm^2")
        
        final_answer_str = f"{a}:{b}:{c}:{area}"
        print(f"\nFinal Answer: <<< {final_answer_str} >>>")
    else:
        # No better solution found with this packing model
        print("\nNo box configuration with a smaller surface area was found.")
        print("The initial 12x12x12 cm box is the most efficient for a simple grid packing.")
        final_answer_str = "0"
        print(f"\nFinal Answer: <<< {final_answer_str} >>>")

# Run the solver
solve_pioneer_probe_packing()