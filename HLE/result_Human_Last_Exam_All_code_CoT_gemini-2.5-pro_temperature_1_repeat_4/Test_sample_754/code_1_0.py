import math

def solve_interstellar_problem():
    # --- Step 0: Define constants ---
    distance_ly = 10.0  # light-years
    speed_c = 0.02     # as a fraction of speed of light 'c'
    ball_radius = 2.0  # cm
    ball_diameter = ball_radius * 2
    initial_energy_per_ball = 100.0  # MJ
    leak_rate_per_year = 0.001  # 0.1%
    container_dims = [12.0, 11.0, 11.0]  # cm
    required_energy = 1000.0  # MJ

    # --- Step 1: Calculate travel time ---
    travel_time_years = distance_ly / speed_c
    print(f"The probe will travel for {travel_time_years} years.")

    # --- Step 2: Determine ball placement and capacity ---
    # For a simple cubic lattice packing
    nx = math.floor((container_dims[0] - ball_diameter) / ball_diameter) + 1
    ny = math.floor((container_dims[1] - ball_diameter) / ball_diameter) + 1
    nz = math.floor((container_dims[2] - ball_diameter) / ball_diameter) + 1
    total_balls = nx * ny * nz
    print(f"The container can hold a maximum of {nx}x{ny}x{nz} = {total_balls} balls.")

    # --- Step 3: Differentiate inner and outer balls ---
    # Center coordinates for the packed balls
    x_centers = [ball_radius + i * ball_diameter for i in range(nx)]
    y_centers = [ball_radius + i * ball_diameter for i in range(ny)]
    z_centers = [ball_radius + i * ball_diameter for i in range(nz)]
    
    num_inner_balls = 0
    for x in x_centers:
        for y in y_centers:
            for z in z_centers:
                # Check if the ball is "inner" (does not touch any wall)
                is_inner = (
                    (x - ball_radius > 0) and (x + ball_radius < container_dims[0]) and
                    (y - ball_radius > 0) and (y + ball_radius < container_dims[1]) and
                    (z - ball_radius > 0) and (z + ball_radius < container_dims[2])
                )
                if is_inner:
                    num_inner_balls += 1
    
    num_outer_balls = total_balls - num_inner_balls
    print(f"Number of inner balls (no leak): {num_inner_balls}")
    print(f"Number of outer balls (leaking): {num_outer_balls}")

    # --- Step 4: Calculate energy upon arrival ---
    # Energy of an inner ball remains the same
    final_energy_inner_ball = initial_energy_per_ball
    
    # Energy of an outer ball decays over time
    final_energy_outer_ball = initial_energy_per_ball * math.pow(1 - leak_rate_per_year, travel_time_years)
    print(f"Energy of one outer ball upon arrival: {final_energy_outer_ball:.2f} MJ")
    
    total_final_energy = (num_inner_balls * final_energy_inner_ball) + (num_outer_balls * final_energy_outer_ball)

    print("\n--- Final Calculation ---")
    print(f"Total Energy = (Inner Balls * Energy) + (Outer Balls * Energy)")
    print(f"Total Energy = ({num_inner_balls} * {final_energy_inner_ball:.2f}) + ({num_outer_balls} * {final_energy_outer_ball:.2f})")
    print(f"Total Energy = {num_inner_balls * final_energy_inner_ball:.2f} MJ + {num_outer_balls * final_energy_outer_ball:.2f} MJ")
    print(f"Total Energy on Arrival: {total_final_energy:.2f} MJ")
    print(f"Required Energy for Operations: {required_energy:.2f} MJ")

    # --- Step 5: Final conclusion ---
    if total_final_energy >= required_energy:
        print("\nThe container is sufficient.")
        final_answer = total_balls
    else:
        print("\nThe energy is NOT sufficient.")
        final_answer = 0
    
    print(f"Final Answer: {final_answer}")
    print(f"<<<{final_answer}>>>")

solve_interstellar_problem()