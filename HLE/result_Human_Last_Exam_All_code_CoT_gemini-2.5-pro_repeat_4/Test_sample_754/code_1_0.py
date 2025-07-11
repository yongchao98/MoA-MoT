import math

def solve_interstellar_problem():
    """
    Solves the Interstellar energy storage problem by simulating ball packing
    and calculating the final energy.
    """
    # --- Problem Constants ---
    DISTANCE_LY = 10.0
    SPEED_C = 0.02
    INITIAL_ENERGY_MJ = 100.0
    LEAK_RATE_PER_YEAR = 0.001
    REQUIRED_ENERGY_MJ = 1000.0
    BALL_RADIUS_CM = 2.0
    BOX_DIMS_CM = [12.0, 11.0, 11.0]
    GRID_STEP_CM = 0.5

    # --- Step 1: Calculate Travel Time ---
    travel_time_years = DISTANCE_LY / SPEED_C
    print(f"1. Calculating travel time...")
    print(f"   - Probe travel time: {DISTANCE_LY} light-years / {SPEED_C}c = {travel_time_years:.0f} years")
    print("-" * 40)

    # --- Step 2: Simulate Ball Packing ---
    print("2. Simulating ball packing to find the maximum number of balls...")
    placed_centers = []
    ball_diameter = 2 * BALL_RADIUS_CM
    # Use squared distance for efficiency
    dist_sq_min = ball_diameter**2

    # Define the valid region for ball centers
    x_min, y_min, z_min = BALL_RADIUS_CM, BALL_RADIUS_CM, BALL_RADIUS_CM
    x_max = BOX_DIMS_CM[0] - BALL_RADIUS_CM
    y_max = BOX_DIMS_CM[1] - BALL_RADIUS_CM
    z_max = BOX_DIMS_CM[2] - BALL_RADIUS_CM

    # Generate all possible grid points for centers
    x_coords = [x_min + i * GRID_STEP_CM for i in range(int(round((x_max - x_min) / GRID_STEP_CM)) + 1)]
    y_coords = [y_min + i * GRID_STEP_CM for i in range(int(round((y_max - y_min) / GRID_STEP_CM)) + 1)]
    z_coords = [z_min + i * GRID_STEP_CM for i in range(int(round((z_max - z_min) / GRID_STEP_CM)) + 1)]
    
    # Greedy packing algorithm (iterating z, then y, then x is a deterministic choice)
    for z in z_coords:
        for y in y_coords:
            for x in x_coords:
                new_center = (x, y, z)
                can_place = True
                for placed_center in placed_centers:
                    # Calculate squared distance to avoid sqrt
                    dist_sq = (
                        (new_center[0] - placed_center[0]) ** 2
                        + (new_center[1] - placed_center[1]) ** 2
                        + (new_center[2] - placed_center[2]) ** 2
                    )
                    # Add a small tolerance for floating point comparisons
                    if dist_sq < dist_sq_min - 1e-9:
                        can_place = False
                        break
                if can_place:
                    placed_centers.append(new_center)

    total_balls = len(placed_centers)
    print(f"   - A greedy packing algorithm fits a maximum of {total_balls} balls.")
    print("-" * 40)

    # --- Step 3: Categorize Balls (Touching vs. Non-touching) ---
    print("3. Categorizing balls based on their position...")
    touching_balls = 0
    for center in placed_centers:
        cx, cy, cz = center
        # A ball is touching if its center is at the boundary of the valid region
        is_touching = (
            math.isclose(cx, x_min) or math.isclose(cx, x_max) or
            math.isclose(cy, y_min) or math.isclose(cy, y_max) or
            math.isclose(cz, z_min) or math.isclose(cz, z_max)
        )
        if is_touching:
            touching_balls += 1
    
    non_touching_balls = total_balls - touching_balls
    print(f"   - Number of balls touching container walls: {touching_balls}")
    print(f"   - Number of balls not touching walls: {non_touching_balls}")
    print("-" * 40)

    # --- Step 4: Calculate Final Energy ---
    print("4. Calculating total energy upon arrival...")
    energy_after_leak = INITIAL_ENERGY_MJ * (1 - LEAK_RATE_PER_YEAR)**travel_time_years
    energy_no_leak = INITIAL_ENERGY_MJ

    total_final_energy = (touching_balls * energy_after_leak) + (non_touching_balls * energy_no_leak)
    
    # Output the components of the final energy equation
    print(f"   - Energy per leaking ball: {INITIAL_ENERGY_MJ:.1f} * (1 - {LEAK_RATE_PER_YEAR}) ^ {travel_time_years:.0f} = {energy_after_leak:.2f} MJ")
    print(f"   - Energy per non-leaking ball: {energy_no_leak:.2f} MJ")
    print("\n   Final Energy Equation:")
    print(f"   Total Energy = ({touching_balls} balls * {energy_after_leak:.2f} MJ) + ({non_touching_balls} balls * {energy_no_leak:.2f} MJ)")
    print(f"   Total Energy = {touching_balls * energy_after_leak:.2f} MJ + {non_touching_balls * energy_no_leak:.2f} MJ = {total_final_energy:.2f} MJ")
    print("-" * 40)

    # --- Step 5: Determine Final Answer ---
    print("5. Final Conclusion:")
    is_sufficient = total_final_energy >= REQUIRED_ENERGY_MJ
    final_answer = total_balls if is_sufficient else 0

    if is_sufficient:
        print(f"   - The total available energy ({total_final_energy:.2f} MJ) IS SUFFICIENT for the required {REQUIRED_ENERGY_MJ:.2f} MJ.")
    else:
        print(f"   - The total available energy ({total_final_energy:.2f} MJ) IS NOT SUFFICIENT for the required {REQUIRED_ENERGY_MJ:.2f} MJ.")
    
    print(f"\n   Is the container enough? {'Yes' if is_sufficient else 'No'}.")
    print(f"   Result (Maximal number of balls if sufficient, else 0): {final_answer}")
    
    return final_answer

# Run the solution and print the final answer in the required format
final_answer = solve_interstellar_problem()
print(f"\n<<<{final_answer}>>>")
