import math

def solve_interstellar_problem():
    """
    Calculates if the energy container is sufficient for the Pioneer probe's mission.
    """
    # Mission parameters
    distance_ly = 10.0
    speed_fraction_c = 0.02
    required_energy_mj = 1000.0

    # Energy ball parameters
    ball_radius_cm = 2.0
    initial_energy_per_ball_mj = 100.0
    leak_rate_per_year = 0.001

    # Container parameters
    box_x, box_y, box_z = 12.0, 11.0, 11.0
    
    # Step 1: Calculate travel time
    travel_time_years = distance_ly / speed_fraction_c

    # Step 2: Calculate the maximum number of balls that can fit
    ball_diameter_cm = 2 * ball_radius_cm
    
    # Calculate how many balls fit along each dimension in a grid layout
    n_x = math.floor((box_x - 2 * ball_radius_cm) / ball_diameter_cm) + 1
    n_y = math.floor((box_y - 2 * ball_radius_cm) / ball_diameter_cm) + 1
    n_z = math.floor((box_z - 2 * ball_radius_cm) / ball_diameter_cm) + 1
    total_balls = n_x * n_y * n_z

    # Step 3: Identify surface vs. inner balls
    # Generate the coordinates of the ball centers
    coords_x = [ball_radius_cm + i * ball_diameter_cm for i in range(n_x)]
    coords_y = [ball_radius_cm + i * ball_diameter_cm for i in range(n_y)]
    coords_z = [ball_radius_cm + i * ball_diameter_cm for i in range(n_z)]

    num_surface_balls = 0
    # Iterate through all ball positions to check if they touch a wall
    for cx in coords_x:
        for cy in coords_y:
            for cz in coords_z:
                is_surface = (
                    (cx - ball_radius_cm < 1e-9) or (cx + ball_radius_cm > box_x - 1e-9) or
                    (cy - ball_radius_cm < 1e-9) or (cy + ball_radius_cm > box_y - 1e-9) or
                    (cz - ball_radius_cm < 1e-9) or (cz + ball_radius_cm > box_z - 1e-9)
                )
                if is_surface:
                    num_surface_balls += 1
    
    num_inner_balls = total_balls - num_surface_balls

    # Step 4: Calculate the total energy upon arrival
    # Energy of a surface ball after leaking
    final_energy_surface_ball = initial_energy_per_ball_mj * ((1 - leak_rate_per_year) ** travel_time_years)
    
    # Energy of an inner ball does not change
    final_energy_inner_ball = initial_energy_per_ball_mj

    # Total energy is the sum from all balls
    total_final_energy = (num_surface_balls * final_energy_surface_ball) + (num_inner_balls * final_energy_inner_ball)

    # Step 5: Output the calculation and the final decision
    print("Final energy calculation:")
    print(f"{num_surface_balls} * {final_energy_surface_ball:.2f} + {num_inner_balls} * {final_energy_inner_ball:.2f} = {total_final_energy:.2f} MJ")
    print(f"Energy required upon arrival: {required_energy_mj:.2f} MJ")

    if total_final_energy >= required_energy_mj:
        print("Conclusion: The stored energy is sufficient.")
        final_answer = total_balls
    else:
        print("Conclusion: The stored energy is NOT sufficient.")
        final_answer = 0
    
    return final_answer

# Run the solver and print the final answer in the required format
final_answer = solve_interstellar_problem()
print(f"<<<{final_answer}>>>")
