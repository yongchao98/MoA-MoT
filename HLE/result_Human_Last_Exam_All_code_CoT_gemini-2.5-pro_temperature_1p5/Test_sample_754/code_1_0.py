import math

def solve_interstellar_problem():
    """
    Solves the Interstellar energy storage problem by calculating if the container
    can hold enough energy for the Pioneer probe's mission.
    """
    # --- Problem Parameters ---
    distance_ly = 10.0
    probe_speed_c = 0.02
    initial_ball_energy = 100.0  # MJ
    leak_rate_per_year = 0.001
    required_total_energy = 1000.0  # MJ
    container_dim = (12, 11, 11)  # cm
    ball_radius = 2.0  # cm

    # --- Step 1: Calculate travel time ---
    print("--- Step 1: Calculate Travel Time ---")
    travel_time_years = distance_ly / probe_speed_c
    print(f"The journey to Pandora will take {distance_ly} light-years / {probe_speed_c}c = {travel_time_years:.0f} years.\n")

    # --- Step 2: Calculate energy remaining in a leaking ball ---
    print("--- Step 2: Calculate Energy in a Leaking Ball ---")
    energy_after_leak = initial_ball_energy * (1 - leak_rate_per_year) ** travel_time_years
    print(f"An energy ball touching the container wall will have {initial_ball_energy:.2f} * (1 - {leak_rate_per_year}) ^ {travel_time_years:.0f} = {energy_after_leak:.2f} MJ of energy upon arrival.")
    print(f"An internal ball (not touching the wall) will retain its full {initial_ball_energy:.2f} MJ.\n")

    # --- Step 3: Calculate the maximum number of balls that fit ---
    print("--- Step 3: Calculate Maximum Ball Capacity ---")
    ball_diameter = 2 * ball_radius
    
    # Using a simple cubic lattice packing
    n_x = math.floor((container_dim[0] - ball_diameter) / ball_diameter) + 1
    n_y = math.floor((container_dim[1] - ball_diameter) / ball_diameter) + 1
    n_z = math.floor((container_dim[2] - ball_diameter) / ball_diameter) + 1
    max_balls = n_x * n_y * n_z
    
    print(f"The container can fit {n_x} balls along the X-axis (12cm).")
    print(f"The container can fit {n_y} balls along the Y-axis (11cm).")
    print(f"The container can fit {n_z} balls along the Z-axis (11cm).")
    print(f"Maximum number of balls that can be stored is {n_x} * {n_y} * {n_z} = {max_balls}.\n")

    # --- Step 4: Calculate the number of internal and touching balls ---
    print("--- Step 4: Count Internal vs. Touching Balls ---")
    
    # Generate the actual coordinates for the ball centers
    x_coords = [ball_radius + i * ball_diameter for i in range(n_x)]
    y_coords = [ball_radius + i * ball_diameter for i in range(n_y)]
    z_coords = [ball_radius + i * ball_diameter for i in range(n_z)]

    # A ball is touching if its center is at a coordinate that places its surface on the container wall
    touching_x_vals = {ball_radius, container_dim[0] - ball_radius}
    touching_y_vals = {ball_radius, container_dim[1] - ball_radius}
    touching_z_vals = {ball_radius, container_dim[2] - ball_radius}
    
    # Count how many coordinate options are internal (not touching)
    num_internal_x = sum(1 for x in x_coords if x not in touching_x_vals)
    num_internal_y = sum(1 for y in y_coords if y not in touching_y_vals)
    num_internal_z = sum(1 for z in z_coords if z not in touching_z_vals)
    
    num_internal_balls = num_internal_x * num_internal_y * num_internal_z
    num_touching_balls = max_balls - num_internal_balls
    
    print(f"A ball is 'internal' if its center is not on the edge of the placement grid.")
    print(f"Number of internal balls = {num_internal_balls}")
    print(f"Number of touching balls = {max_balls} - {num_internal_balls} = {num_touching_balls}\n")

    # --- Step 5: Calculate the total available energy ---
    print("--- Step 5: Calculate Total Available Energy ---")
    total_energy_available = (num_internal_balls * initial_ball_energy) + (num_touching_balls * energy_after_leak)
    
    print("The final equation for total energy is:")
    print(f"Total Energy = (Number of Internal Balls * Energy per Internal Ball) + (Number of Touching Balls * Energy per Touching Ball)")
    print(f"Total Energy = ({num_internal_balls} * {initial_ball_energy:.2f}) + ({num_touching_balls} * {energy_after_leak:.2f})")
    print(f"Total Energy = {num_internal_balls * initial_ball_energy:.2f} MJ + {num_touching_balls * energy_after_leak:.2f} MJ = {total_energy_available:.2f} MJ\n")

    # --- Step 6: Final Conclusion ---
    print("--- Step 6: Final Conclusion ---")
    print(f"The energy required for operations is {required_total_energy:.2f} MJ.")
    print(f"The total energy available from {max_balls} balls is {total_energy_available:.2f} MJ.")

    if total_energy_available >= required_total_energy:
        print("This is SUFFICIENT. The probe can operate.")
        final_answer = max_balls
    else:
        print("This is NOT SUFFICIENT. The probe cannot operate.")
        final_answer = 0
        
    print(f"\nThe answer is the maximal number of balls if the energy is sufficient, otherwise 0.")
    print(f"<<<{final_answer}>>>")

# Run the solution
solve_interstellar_problem()