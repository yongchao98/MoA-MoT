import math

def solve_interstellar_problem():
    """
    Calculates if the Pioneer probe's energy container is sufficient for its mission.
    """

    # --- Step 1: Define constants and calculate travel time ---
    distance_ly = 10.0  # light-years
    speed_c = 0.02  # as a fraction of the speed of light
    travel_time_years = distance_ly / speed_c

    # --- Step 2: Define physical parameters and calculate ball capacity ---
    container_dims = (12.0, 11.0, 11.0)  # in cm
    ball_radius = 2.0  # in cm
    ball_diameter = ball_radius * 2
    
    # Max balls in a simple grid packing
    nx = math.floor(container_dims[0] / ball_diameter)
    ny = math.floor(container_dims[1] / ball_diameter)
    nz = math.floor(container_dims[2] / ball_diameter)
    max_balls = nx * ny * nz

    print(f"--- Mission Analysis ---")
    print(f"1. Travel Time Calculation")
    print(f"   The probe travels 10 light-years at 0.02c.")
    print(f"   Time = Distance / Speed = {distance_ly} ly / {speed_c}c = {int(travel_time_years)} years\n")
    
    print(f"2. Container Capacity Calculation")
    print(f"   Container (cm): {container_dims[0]} x {container_dims[1]} x {container_dims[2]}")
    print(f"   Energy ball diameter: {ball_diameter} cm")
    print(f"   Balls per dimension (x, y, z): {nx}, {ny}, {nz}")
    print(f"   Maximum number of balls: {nx} * {ny} * {nz} = {max_balls}\n")

    # --- Step 3: Classify inner and outer balls ---
    num_inner_balls = 0
    ball_centers = []
    # Generate centers for all packed balls
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                # Place first ball center at (radius, radius, radius)
                # and subsequent balls at intervals of diameter
                cx = ball_radius + i * ball_diameter
                cy = ball_radius + j * ball_diameter
                cz = ball_radius + k * ball_diameter
                
                # Check if the ball touches any container wall
                is_outer = (
                    (cx - ball_radius <= 0) or (cx + ball_radius >= container_dims[0]) or
                    (cy - ball_radius <= 0) or (cy + ball_radius >= container_dims[1]) or
                    (cz - ball_radius <= 0) or (cz + ball_radius >= container_dims[2])
                )
                
                if not is_outer:
                    num_inner_balls += 1
    
    num_outer_balls = max_balls - num_inner_balls
    
    print(f"3. Inner and Outer Ball Classification")
    print(f"   An 'outer' ball is one touching any surface of the container.")
    print(f"   Based on the packing configuration:")
    print(f"   Number of outer (leaking) balls: {num_outer_balls}")
    print(f"   Number of inner (stable) balls: {num_inner_balls}\n")

    # --- Step 4 & 5: Calculate total energy ---
    initial_energy_mj = 100.0
    leak_rate_per_year = 0.001  # 0.1%
    required_energy_mj = 1000.0

    energy_per_leaking_ball = initial_energy_mj * ((1 - leak_rate_per_year) ** travel_time_years)
    
    total_energy_from_outer = num_outer_balls * energy_per_leaking_ball
    total_energy_from_inner = num_inner_balls * initial_energy_mj
    total_final_energy = total_energy_from_inner + total_energy_from_outer
    
    print(f"4. Energy Calculation upon Arrival")
    print(f"   Initial energy per ball: {initial_energy_mj} MJ")
    print(f"   Energy of 1 outer ball after {int(travel_time_years)} years = {initial_energy_mj} * (1 - {leak_rate_per_year}) ^ {int(travel_time_years)} = {energy_per_leaking_ball:.2f} MJ")
    print(f"   Total energy = ({num_inner_balls} inner * {initial_energy_mj} MJ) + ({num_outer_balls} outer * {energy_per_leaking_ball:.2f} MJ)")
    print(f"   Total available energy = {total_energy_from_inner:.2f} MJ + {total_energy_from_outer:.2f} MJ = {total_final_energy:.2f} MJ\n")

    # --- Step 6: Compare and conclude ---
    is_sufficient = total_final_energy >= required_energy_mj
    final_answer = max_balls if is_sufficient else 0

    print(f"5. Final Conclusion")
    print(f"   Required energy for operations: {required_energy_mj} MJ")
    print(f"   Is {total_final_energy:.2f} MJ >= {required_energy_mj} MJ? {'Yes' if is_sufficient else 'No'}")
    
    if is_sufficient:
        print(f"\nThe energy is sufficient. The answer is the maximal number of balls.")
    else:
        print(f"\nThe available energy is not sufficient for the mission.")
    
    # --- Final formatted answer ---
    print(f"\n<<<>>>") # Required format start
    print(final_answer)
    print(f"<<<>>>") # Required format end

solve_interstellar_problem()