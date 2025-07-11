import math

def solve_interstellar_problem():
    """
    Calculates if the energy container is sufficient for the Pioneer probe's mission.
    """
    # --- Constants ---
    distance_ly = 10.0
    speed_c = 0.02
    initial_energy_mj = 100.0
    leakage_rate_yearly = 0.001
    ball_radius_cm = 2.0
    box_dims_cm = (12.0, 11.0, 11.0)
    required_energy_mj = 1000.0
    
    # --- Step 1: Calculate travel time ---
    travel_time_years = distance_ly / speed_c

    # --- Step 2: Calculate energy of a leaking ball ---
    energy_leaking_ball = initial_energy_mj * (1 - leakage_rate_yearly) ** travel_time_years
    energy_internal_ball = initial_energy_mj

    # --- Step 3: Determine maximum number of balls (simple grid packing) ---
    ball_diameter = 2 * ball_radius_cm
    max_balls = [0, 0, 0]
    for i in range(3):
        dim_length = box_dims_cm[i]
        # The valid range for centers is from radius to dim_length - radius
        valid_range_length = dim_length - 2 * ball_radius_cm
        if valid_range_length >= 0:
            # Number of items you can fit in a line is floor(space / item_size) + 1
            max_balls[i] = math.floor(valid_range_length / ball_diameter) + 1
    
    total_balls_packed = max_balls[0] * max_balls[1] * max_balls[2]

    # --- Step 4: Count touching and non-touching balls ---
    # We generate the center coordinates for the packed balls
    # We start placing balls at the minimum coordinate, which is 'radius'
    x_centers = [ball_radius_cm + i * ball_diameter for i in range(max_balls[0])]
    y_centers = [ball_radius_cm + i * ball_diameter for i in range(max_balls[1])]
    z_centers = [ball_radius_cm + i * ball_diameter for i in range(max_balls[2])]
    
    num_touching = 0
    num_internal = 0

    for x in x_centers:
        for y in y_centers:
            for z in z_centers:
                # A ball is touching if its surface touches any of the 6 container walls
                is_touching = ( (x - ball_radius_cm <= 0) or (x + ball_radius_cm >= box_dims_cm[0]) or
                                (y - ball_radius_cm <= 0) or (y + ball_radius_cm >= box_dims_cm[1]) or
                                (z - ball_radius_cm <= 0) or (z + ball_radius_cm >= box_dims_cm[2]) )
                
                if is_touching:
                    num_touching += 1
                else:
                    num_internal += 1

    # --- Step 5: Calculate total available energy ---
    total_energy = (num_internal * energy_internal_ball) + (num_touching * energy_leaking_ball)

    # --- Print the process and the equation ---
    print(f"Maximum number of balls that can be packed: {total_balls_packed}")
    print(f"Number of internal (non-leaking) balls: {num_internal}")
    print(f"Number of touching (leaking) balls: {num_touching}")
    print(f"Energy per leaking ball after {int(travel_time_years)} years: {energy_leaking_ball:.2f} MJ")
    print("\nFinal energy equation:")
    print(f"Total Energy = ({num_internal} * {initial_energy_mj:.0f} MJ) + ({num_touching} * {energy_leaking_ball:.2f} MJ)")
    print(f"Total Energy = {num_internal * initial_energy_mj:.2f} MJ + {num_touching * energy_leaking_ball:.2f} MJ = {total_energy:.2f} MJ")
    
    # --- Step 6: Final Decision ---
    if total_energy >= required_energy_mj:
        final_answer = total_balls_packed
    else:
        final_answer = 0
        
    print(f"\nTotal available energy ({total_energy:.2f} MJ) is less than the required energy ({required_energy_mj:.0f} MJ).")
    print(f"Conclusion: The energy container is not sufficient.")

    # --- Output the final answer ---
    print(f"\nFinal Answer: {final_answer}")
    return final_answer

# Run the simulation and print the final answer in the required format
final_result = solve_interstellar_problem()
print(f"<<<{final_result}>>>")