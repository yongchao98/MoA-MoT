import numpy as np

def solve_interstellar_problem():
    """
    Solves the Interstellar energy storage problem.
    """
    # --- Part 1: Define constants and calculate travel time ---
    distance_ly = 10.0
    speed_c = 0.02
    travel_time = distance_ly / speed_c

    box_dims = (12.0, 11.0, 11.0)
    radius = 2.0
    grid_step = 0.5
    diameter = 2 * radius
    diameter_sq = diameter**2
    
    initial_energy_per_ball = 100.0  # MJ
    leak_rate = 0.001  # 0.1% per year
    required_energy = 1000.0  # MJ
    
    print("Step 1: Calculating the travel time.")
    print(f"The probe travels 10 light-years at 0.02c. Time = Distance / Speed.")
    print(f"Travel Time = {distance_ly} light-years / {speed_c}c = {int(travel_time)} years.\n")

    # --- Part 2: Greedy Packing Algorithm to find max balls ---
    print("Step 2: Determining the maximum number of energy balls that can be packed.")
    print("We use a greedy algorithm to place balls (radius 2cm) in the 12x11x11 cm box.")
    print("The centers must be on a 0.5cm grid and not overlap.\n")

    placed_centers = []
    
    x_min, y_min, z_min = radius, radius, radius
    x_max = box_dims[0] - radius
    y_max = box_dims[1] - radius
    z_max = box_dims[2] - radius
    
    x_coords = np.arange(x_min, x_max + grid_step * 0.5, grid_step)
    y_coords = np.arange(y_min, y_max + grid_step * 0.5, grid_step)
    z_coords = np.arange(z_min, z_max + grid_step * 0.5, grid_step)

    epsilon = 1e-9
    for z in z_coords:
        for y in y_coords:
            for x in x_coords:
                candidate_center = (x, y, z)
                is_valid = True
                for placed_center in placed_centers:
                    dist_sq = sum([(c1 - c2)**2 for c1, c2 in zip(candidate_center, placed_center)])
                    if dist_sq < diameter_sq - epsilon:
                        is_valid = False
                        break
                if is_valid:
                    placed_centers.append(candidate_center)

    num_total_balls = len(placed_centers)
    print(f"The maximum number of balls that can be packed is: {num_total_balls}\n")

    # --- Part 3: Calculate energy upon arrival ---
    print("Step 3: Calculating the total energy upon arrival.")
    
    num_leaking_balls = 0
    fuzz = 1e-6
    leaking_boundaries = {
        'x': [x_min, x_max],
        'y': [y_min, y_max],
        'z': [z_min, z_max]
    }
    
    for center in placed_centers:
        cx, cy, cz = center
        if (abs(cx - leaking_boundaries['x'][0]) < fuzz or abs(cx - leaking_boundaries['x'][1]) < fuzz or
            abs(cy - leaking_boundaries['y'][0]) < fuzz or abs(cy - leaking_boundaries['y'][1]) < fuzz or
            abs(cz - leaking_boundaries['z'][0]) < fuzz or abs(cz - leaking_boundaries['z'][1]) < fuzz):
            num_leaking_balls += 1
            
    num_inner_balls = num_total_balls - num_leaking_balls
    
    print(f"Out of {num_total_balls} balls, {num_leaking_balls} are touching the container walls (leaking) and {num_inner_balls} are not (inner).")

    energy_leaking_ball = initial_energy_per_ball * ((1 - leak_rate) ** travel_time)
    print(f"Energy of one leaking ball after {int(travel_time)} years = 100 MJ * (1 - {leak_rate})^{int(travel_time)} = {energy_leaking_ball:.2f} MJ.")
    
    total_energy = (num_inner_balls * initial_energy_per_ball) + (num_leaking_balls * energy_leaking_ball)
    
    print("\nFinal Energy Calculation:")
    print(f"Total Energy = ({num_inner_balls} inner balls * {initial_energy_per_ball:.2f} MJ/ball) + ({num_leaking_balls} leaking balls * {energy_leaking_ball:.2f} MJ/ball)")
    print(f"Total Energy = {num_inner_balls * initial_energy_per_ball:.2f} MJ + {num_leaking_balls * energy_leaking_ball:.2f} MJ = {total_energy:.2f} MJ\n")
    
    # --- Part 4: Final Conclusion ---
    print("Step 4: Conclusion.")
    if total_energy >= required_energy:
        print(f"The calculated total energy ({total_energy:.2f} MJ) is greater than or equal to the required {required_energy:.2f} MJ.")
        print("Therefore, the energy container is sufficient.")
        final_answer = num_total_balls
    else:
        print(f"The calculated total energy ({total_energy:.2f} MJ) is less than the required {required_energy:.2f} MJ.")
        print("Therefore, the energy container is NOT sufficient.")
        final_answer = 0
        
    return final_answer

if __name__ == '__main__':
    answer = solve_interstellar_problem()
    print(f"\n<<<{answer}>>>")