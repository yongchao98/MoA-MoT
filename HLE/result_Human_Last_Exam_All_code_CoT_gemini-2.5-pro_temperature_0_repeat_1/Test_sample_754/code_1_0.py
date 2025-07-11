import math

def solve_interstellar_problem():
    """
    Calculates if the energy container is sufficient for the Pioneer probe's mission to Pandora.
    """
    # Step 1: Define all the constants from the problem description.
    distance_ly = 10.0  # light-years
    speed_c = 0.02  # fraction of the speed of light
    initial_energy_per_ball_MJ = 100.0
    leak_rate_per_year = 0.001  # 0.1%
    required_energy_MJ = 1000.0
    ball_radius_cm = 2.0
    ball_diameter_cm = 2 * ball_radius_cm
    container_dims_cm = (12.0, 11.0, 11.0)
    L, W, H = container_dims_cm

    # Step 2: Calculate the travel time in years.
    travel_time_years = distance_ly / speed_c

    # Step 3: Calculate the maximum number of balls that can fit in a grid layout.
    # The number of balls along one dimension is floor((Length - 2*radius) / diameter) + 1.
    nx = math.floor((L - 2 * ball_radius_cm) / ball_diameter_cm) + 1
    ny = math.floor((W - 2 * ball_radius_cm) / ball_diameter_cm) + 1
    nz = math.floor((H - 2 * ball_radius_cm) / ball_diameter_cm) + 1
    max_balls = nx * ny * nz

    # Step 4: Categorize balls into "touching" and "internal".
    # We iterate through each ball's position in the grid to check if it touches a wall.
    num_internal_balls = 0
    num_touching_balls = 0
    
    # Define the coordinates of the centers of the balls in the grid.
    # The first ball's center is at `radius`, and subsequent centers are `diameter` apart.
    centers_x = [ball_radius_cm + i * ball_diameter_cm for i in range(nx)]
    centers_y = [ball_radius_cm + i * ball_diameter_cm for i in range(ny)]
    centers_z = [ball_radius_cm + i * ball_diameter_cm for i in range(nz)]

    for cx in centers_x:
        for cy in centers_y:
            for cz in centers_z:
                # A ball is touching if its surface is at or beyond the container's boundary.
                is_touching = (
                    (cx - ball_radius_cm <= 0) or
                    (cx + ball_radius_cm >= L) or
                    (cy - ball_radius_cm <= 0) or
                    (cy + ball_radius_cm >= W) or
                    (cz - ball_radius_cm <= 0) or
                    (cz + ball_radius_cm >= H)
                )
                if is_touching:
                    num_touching_balls += 1
                else:
                    num_internal_balls += 1

    # Step 5: Calculate the total energy upon arrival.
    # Energy of a single touching ball after decay.
    energy_per_touching_ball = initial_energy_per_ball_MJ * (1 - leak_rate_per_year) ** travel_time_years
    # Internal balls do not leak energy.
    energy_per_internal_ball = initial_energy_per_ball_MJ
    
    total_final_energy_MJ = (num_touching_balls * energy_per_touching_ball) + (num_internal_balls * energy_per_internal_ball)

    # Step 6: Compare with requirement and determine the final answer.
    is_sufficient = total_final_energy_MJ >= required_energy_MJ
    final_answer = max_balls if is_sufficient else 0

    # --- Output ---
    print(f"The probe will travel for {int(travel_time_years)} years.")
    print(f"The {L}x{W}x{H} cm container can hold a maximum of {max_balls} balls ({nx}x{ny}x{nz}).")
    print(f"Of these, {num_touching_balls} balls touch the container walls and will leak energy.")
    print(f"There are {num_internal_balls} internal balls that will not leak energy.")
    print("\n--- Energy Calculation ---")
    
    # Print the final equation with all its numbers, as requested.
    print("Final Energy Equation:")
    print(f"{num_touching_balls} * {initial_energy_per_ball_MJ} * (1 - {leak_rate_per_year})^{int(travel_time_years)} + {num_internal_balls} * {initial_energy_per_ball_MJ} = {total_final_energy_MJ:.2f} MJ")
    
    print(f"\nTotal energy on arrival: {total_final_energy_MJ:.2f} MJ")
    print(f"Required energy for operations: {required_energy_MJ:.2f} MJ")

    if is_sufficient:
        print("\nResult: The available energy is SUFFICIENT for the mission.")
    else:
        print("\nResult: The available energy is NOT SUFFICIENT for the mission.")
    
    # The final answer in the required format.
    print(f"\nIs the container enough? (Max balls if Yes, 0 if No)")
    print(f"<<<{final_answer}>>>")

# Execute the function
solve_interstellar_problem()