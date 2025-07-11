import math

def solve_interstellar_problem():
    """
    Solves the Interstellar energy storage problem by calculating if the container
    can hold enough energy balls for the Pioneer probe's mission.
    """
    # --- Step 1: Define constants and calculate travel time ---
    distance_ly = 10.0
    speed_frac_c = 0.02
    initial_energy_per_ball_mj = 100.0
    leak_rate_per_year = 0.001
    required_energy_mj = 1000.0

    box_dims = {'x': 12.0, 'y': 11.0, 'z': 11.0}
    ball_radius_cm = 2.0
    ball_diameter_cm = 2.0 * ball_radius_cm

    travel_time_years = distance_ly / speed_frac_c

    print(f"Pioneer's travel time to Pandora:")
    print(f"Equation: {distance_ly} light-years / {speed_frac_c}c = {travel_time_years} years\n")

    # --- Step 2: Calculate maximum ball capacity ---
    # The first ball's center is at `ball_radius_cm` from the wall.
    # The number of balls along a dimension is floor(space / spacing) + 1.
    balls_in_dim = {}
    for dim in ['x', 'y', 'z']:
        available_space = box_dims[dim] - 2 * ball_radius_cm
        balls_in_dim[dim] = math.floor(available_space / ball_diameter_cm) + 1
    
    total_balls = balls_in_dim['x'] * balls_in_dim['y'] * balls_in_dim['z']

    print(f"Maximum number of energy balls in the {box_dims['x']}x{box_dims['y']}x{box_dims['z']} cm container:")
    print(f"Balls along x-axis = {balls_in_dim['x']}")
    print(f"Balls along y-axis = {balls_in_dim['y']}")
    print(f"Balls along z-axis = {balls_in_dim['z']}")
    print(f"Total number of balls = {balls_in_dim['x']} * {balls_in_dim['y']} * {balls_in_dim['z']} = {total_balls}\n")

    # --- Step 3: Identify touching and non-touching balls ---
    # Generate the center coordinates for balls in each dimension
    centers = {}
    for dim in ['x', 'y', 'z']:
        centers[dim] = [ball_radius_cm + i * ball_diameter_cm for i in range(balls_in_dim[dim])]

    # Count how many center positions are "non-touching" in each dimension
    non_touching_counts = {}
    for dim in ['x', 'y', 'z']:
        count = 0
        for center_pos in centers[dim]:
            # A position is non-touching if it's clear of both walls
            if (center_pos - ball_radius_cm > 0) and (center_pos + ball_radius_cm < box_dims[dim]):
                count += 1
        non_touching_counts[dim] = count
        
    non_touching_balls = non_touching_counts['x'] * non_touching_counts['y'] * non_touching_counts['z']
    touching_balls = total_balls - non_touching_balls

    print(f"Counting the balls based on their position:")
    print(f"Number of non-touching balls (inner core) = {non_touching_counts['x']} * {non_touching_counts['y']} * {non_touching_counts['z']} = {non_touching_balls}")
    print(f"Number of touching balls (outer layer) = {total_balls} - {non_touching_balls} = {touching_balls}\n")

    # --- Step 4: Calculate total energy on arrival ---
    energy_per_non_touching_ball = initial_energy_per_ball_mj
    
    print(f"Calculating final energy after {travel_time_years} years:")
    print("Energy of one non-touching ball = 100 MJ")
    print(f"Energy of one touching ball (leaking {leak_rate_per_year*100}% per year):")
    final_energy_per_touching_ball = initial_energy_per_ball_mj * (1 - leak_rate_per_year)**travel_time_years
    print(f"Equation: {initial_energy_per_ball_mj} * (1 - {leak_rate_per_year}) ^ {travel_time_years} = {final_energy_per_touching_ball:.4f} MJ\n")

    total_final_energy = (touching_balls * final_energy_per_touching_ball) + \
                         (non_touching_balls * energy_per_non_touching_ball)
    
    print("Calculating total available energy upon arrival:")
    print(f"Equation: ({touching_balls} touching balls * {final_energy_per_touching_ball:.4f} MJ) + ({non_touching_balls} non-touching balls * {energy_per_non_touching_ball} MJ)")
    print(f"Total energy = {total_final_energy:.4f} MJ\n")
    
    # --- Step 5: Final verdict ---
    print("Final Verdict:")
    print(f"The required energy for operations is {required_energy_mj} MJ.")
    
    if total_final_energy >= required_energy_mj:
        final_answer = total_balls
        print(f"The available energy ({total_final_energy:.2f} MJ) is sufficient.")
    else:
        final_answer = 0
        print(f"The available energy ({total_final_energy:.2f} MJ) is NOT sufficient.")
    
    return final_answer

if __name__ == '__main__':
    result = solve_interstellar_problem()
    print(f"\n<<<Result>>>\n{result}")
    print(f"\n<<<{result}>>>")