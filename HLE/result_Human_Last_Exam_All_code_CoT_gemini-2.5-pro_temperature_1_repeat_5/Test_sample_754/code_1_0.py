import math

def solve_interstellar_problem():
    """
    Calculates if the Pioneer probe has enough energy for its mission to Pandora.
    """
    # Step 1: Define constants from the problem description
    distance_ly = 10.0
    probe_speed_c = 0.02
    required_energy_MJ = 1000.0
    box_dims = (12.0, 11.0, 11.0)  # (L, W, H) in cm
    ball_radius_cm = 2.0
    ball_initial_energy_MJ = 100.0
    energy_leak_rate = 0.001  # 0.1% per year

    # Step 2: Calculate the travel time
    travel_time_years = distance_ly / probe_speed_c
    print(f"Step 1: Calculate travel time")
    print(f"The journey to Pandora will take {distance_ly} light-years / {probe_speed_c}c = {travel_time_years:.0f} years.\n")

    # Step 3: Calculate the maximum number of balls that can be packed
    ball_diameter_cm = 2 * ball_radius_cm
    nx = math.floor(box_dims[0] / ball_diameter_cm)
    ny = math.floor(box_dims[1] / ball_diameter_cm)
    nz = math.floor(box_dims[2] / ball_diameter_cm)
    max_balls = nx * ny * nz
    
    print(f"Step 2: Calculate the maximum number of energy balls in the container")
    print(f"Container dimensions: {box_dims[0]}x{box_dims[1]}x{box_dims[2]} cm. Ball diameter: {ball_diameter_cm} cm.")
    print(f"Balls fit along length (12 cm): floor(12 / 4) = {nx}")
    print(f"Balls fit along width (11 cm): floor(11 / 4) = {ny}")
    print(f"Balls fit along height (11 cm): floor(11 / 4) = {nz}")
    print(f"Maximum number of balls = {nx} * {ny} * {nz} = {max_balls}.\n")

    # Step 4: Determine the number of leaking (touching) and non-leaking balls
    # Balls touch the wall only if the packing is perfectly snug along that dimension.
    # Gap along x-axis = 12 - (3 * 4) = 0. So balls touch the x=0 and x=12 walls.
    # Gap along y/z-axis = 11 - (2 * 4) = 3. So balls do not touch y/z walls.
    # The balls touching the walls are the two layers at the front and back of the x-axis.
    touching_balls = 2 * ny * nz
    non_touching_balls = max_balls - touching_balls
    
    print(f"Step 3: Determine which balls are leaking energy")
    print(f"The packing is only tight along the 12cm dimension, so only balls on the front and back layers touch the container.")
    print(f"Number of touching (leaking) balls = 2 layers * {ny} balls/layer * {nz} balls/layer = {touching_balls}")
    print(f"Number of non-touching (safe) balls = {max_balls} - {touching_balls} = {non_touching_balls}.\n")

    # Step 5: Calculate the total energy upon arrival
    energy_non_leaking_ball = ball_initial_energy_MJ
    energy_leaking_ball = ball_initial_energy_MJ * (1 - energy_leak_rate) ** travel_time_years
    total_final_energy = (touching_balls * energy_leaking_ball) + (non_touching_balls * energy_non_leaking_ball)
    
    print(f"Step 4: Calculate the total available energy upon arrival")
    print(f"Energy of a single non-leaking ball = {energy_non_leaking_ball:.2f} MJ")
    print(f"Energy of a single leaking ball after {travel_time_years:.0f} years = {ball_initial_energy_MJ} * (1 - {energy_leak_rate})^{travel_time_years:.0f} = {energy_leaking_ball:.2f} MJ")
    print("\nFinal total energy equation:")
    print(f"Total Energy = ({touching_balls} leaking balls × {energy_leaking_ball:.2f} MJ) + ({non_touching_balls} safe balls × {energy_non_leaking_ball:.2f} MJ)")
    print(f"Total Energy = {touching_balls * energy_leaking_ball:.2f} MJ + {non_touching_balls * energy_non_leaking_ball:.2f} MJ")
    print(f"Total Energy = {total_final_energy:.2f} MJ.\n")

    # Step 6: Make the final decision
    print(f"Step 5: Compare available energy with required energy")
    print(f"Available energy ({total_final_energy:.2f} MJ) vs. Required energy ({required_energy_MJ:.2f} MJ)")
    
    if total_final_energy >= required_energy_MJ:
        result = max_balls
        print("Conclusion: The energy is sufficient for the mission.")
    else:
        result = 0
        print("Conclusion: The energy is NOT sufficient for the mission.")
    
    print(f"\nFinal Answer: {result}")
    print(f"<<<{result}>>>")

solve_interstellar_problem()