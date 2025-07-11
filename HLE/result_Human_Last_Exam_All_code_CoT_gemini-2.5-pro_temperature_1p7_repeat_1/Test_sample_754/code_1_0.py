import numpy as np

def solve_interstellar_problem():
    """
    Solves the Pioneer probe energy storage problem.
    """
    # --- Input Parameters ---
    distance_ly = 10.0
    speed_fraction_c = 0.02
    initial_energy_mj = 100.0
    leak_rate_per_year = 0.001
    required_energy_mj = 1000.0
    box_dims = np.array([12.0, 11.0, 11.0])
    ball_radius = 2.0
    grid_step = 0.5

    # --- Step 1: Calculate Travel Time ---
    travel_time_years = distance_ly / speed_fraction_c
    print(f"1. Travel Time Calculation:")
    print(f"   Travel time = {distance_ly} light-years / {speed_fraction_c}c = {int(travel_time_years)} years\n")

    # --- Step 2: Calculate Energy Decay for a Single Touching Ball ---
    final_energy_leaking_ball = initial_energy_mj * (1 - leak_rate_per_year) ** travel_time_years
    print(f"2. Energy Decay Calculation (for one touching ball):")
    print(f"   Energy after {int(travel_time_years)} years = {initial_energy_mj} MJ * (1 - {leak_rate_per_year})^{int(travel_time_years)} = {final_energy_leaking_ball:.2f} MJ\n")

    # --- Step 3: Determine Maximum Number of Balls (Sphere Packing) ---
    print("3. Sphere Packing Simulation:")
    
    # Define the valid region for ball centers
    min_coords = ball_radius
    max_coords = box_dims - ball_radius
    ball_diameter = 2 * ball_radius
    min_dist_sq = ball_diameter**2

    # Generate all potential center points on the grid
    x_coords = np.arange(min_coords, max_coords[0] + grid_step, grid_step)
    y_coords = np.arange(min_coords, max_coords[1] + grid_step, grid_step)
    z_coords = np.arange(min_coords, max_coords[2] + grid_step, grid_step)

    potential_centers = []
    for x in x_coords:
        for y in y_coords:
            for z in z_coords:
                potential_centers.append(np.array([x, y, z]))

    # Greedy packing algorithm
    placed_centers = []
    for p_potential in potential_centers:
        is_valid_placement = True
        for p_placed in placed_centers:
            # Using squared distance to avoid sqrt calculation
            dist_sq = np.sum((p_potential - p_placed)**2)
            if dist_sq < min_dist_sq:
                is_valid_placement = False
                break
        if is_valid_placement:
            placed_centers.append(p_potential)
    
    max_balls = len(placed_centers)
    print(f"   Maximal number of balls that can be packed: {max_balls}\n")

    # --- Step 4: Calculate Total Available Energy ---
    print("4. Total Energy Calculation:")
    num_touching = 0
    num_internal = 0
    
    if max_balls > 0:
        for center in placed_centers:
            # A ball is "touching" if its center is on the boundary of the allowed volume
            if (np.isclose(center[0], min_coords) or np.isclose(center[0], max_coords[0]) or
                np.isclose(center[1], min_coords) or np.isclose(center[1], max_coords[1]) or
                np.isclose(center[2], min_coords) or np.isclose(center[2], max_coords[2])):
                num_touching += 1
            else:
                num_internal += 1
    
    print(f"   Number of balls touching the container: {num_touching}")
    print(f"   Number of internal balls (not touching): {num_internal}\n")

    total_final_energy = (num_touching * final_energy_leaking_ball) + (num_internal * initial_energy_mj)
    
    print("   Final equation for total energy:")
    print(f"   Total Energy = {num_touching} * {final_energy_leaking_ball:.2f} MJ + {num_internal} * {initial_energy_mj:.2f} MJ = {total_final_energy:.2f} MJ\n")

    # --- Step 5: Final Conclusion ---
    print("5. Conclusion:")
    is_sufficient = total_final_energy >= required_energy_mj
    if is_sufficient:
        print(f"   {total_final_energy:.2f} MJ is >= {required_energy_mj} MJ.")
        print("   The energy is sufficient for the mission.")
        final_answer = max_balls
    else:
        print(f"   {total_final_energy:.2f} MJ is < {required_energy_mj} MJ.")
        print("   The energy is NOT sufficient for the mission.")
        final_answer = 0
        
    return final_answer

if __name__ == '__main__':
    # Run the simulation and print the final answer in the required format
    answer = solve_interstellar_problem()
    print(f"\nIs the container enough? If Yes, maximal number of balls, otherwise 0.")
    print(f"<<<{answer}>>>")
