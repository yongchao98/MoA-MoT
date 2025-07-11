import math
import numpy as np

def solve():
    # Step 1: Define constants and problem parameters
    distance_ly = 10
    speed_c = 0.02
    travel_time_years = distance_ly / speed_c

    initial_energy_per_ball_MJ = 100
    required_energy_MJ = 1000
    leak_rate_per_year = 0.001

    box_dims = np.array([12, 11, 11])
    ball_radius = 2.0
    min_center_dist = 2 * ball_radius
    coord_step = 0.5

    # Step 2: Determine the valid region for ball centers
    center_min = np.array([ball_radius, ball_radius, ball_radius])
    center_max = box_dims - ball_radius

    # Step 3: Implement the greedy packing algorithm (center-out)
    
    # Generate all possible center points
    potential_centers = []
    x_coords = np.arange(center_min[0], center_max[0] + coord_step, coord_step)
    y_coords = np.arange(center_min[1], center_max[1] + coord_step, coord_step)
    z_coords = np.arange(center_min[2], center_max[2] + coord_step, coord_step)

    for z in z_coords:
        for y in y_coords:
            for x in x_coords:
                potential_centers.append(np.array([x, y, z]))

    # Sort points by distance from the geometric center of the allowed volume
    volume_center = (center_min + center_max) / 2.0
    potential_centers.sort(key=lambda p: np.sum((p - volume_center)**2))

    # Greedily place balls
    placed_centers = []
    min_dist_sq = min_center_dist**2
    for p_center in potential_centers:
        can_place = True
        for placed in placed_centers:
            dist_sq = np.sum((p_center - placed)**2)
            if dist_sq < min_dist_sq:
                can_place = False
                break
        if can_place:
            placed_centers.append(p_center)

    n_total = len(placed_centers)
    print(f"The packing algorithm placed a maximum of {n_total} balls.")

    # Step 4: Count touching vs. non-touching balls
    n_touching = 0
    for center in placed_centers:
        # A ball is touching if its center is at the boundary of the allowed center region
        if (center[0] == center_min[0] or center[0] == center_max[0] or
            center[1] == center_min[1] or center[1] == center_max[1] or
            center[2] == center_min[2] or center[2] == center_max[2]):
            n_touching += 1
    
    n_non_touching = n_total - n_touching
    print(f"Number of balls touching the container: {n_touching}")
    print(f"Number of balls not touching (internal): {n_non_touching}")

    # Step 5: Calculate the final energy
    energy_factor = (1 - leak_rate_per_year) ** travel_time_years
    energy_per_touching_ball = initial_energy_per_ball_MJ * energy_factor
    energy_per_non_touching_ball = initial_energy_per_ball_MJ

    total_final_energy = (n_non_touching * energy_per_non_touching_ball) + \
                         (n_touching * energy_per_touching_ball)

    print(f"\nTravel time: {int(travel_time_years)} years")
    print(f"Energy of a non-touching ball upon arrival: {energy_per_non_touching_ball:.2f} MJ")
    print(f"Energy of a touching ball upon arrival: {energy_per_touching_ball:.2f} MJ")
    
    print("\nFinal energy calculation:")
    print(f"Total Energy = ({n_non_touching} non-touching balls * {energy_per_non_touching_ball:.2f} MJ) + ({n_touching} touching balls * {energy_per_touching_ball:.2f} MJ)")
    print(f"Total Energy = {n_non_touching * energy_per_non_touching_ball:.2f} MJ + {n_touching * energy_per_touching_ball:.2f} MJ = {total_final_energy:.2f} MJ")

    # Step 6: Final decision
    print(f"\nRequired energy for operations: {required_energy_MJ} MJ")
    if total_final_energy >= required_energy_MJ:
        print("The total energy is sufficient.")
        final_answer = n_total
    else:
        print("The total energy is NOT sufficient.")
        final_answer = 0
        
    print(f"\nFinal Answer: {final_answer}")
    return final_answer

# Execute the solution
final_answer = solve()
print(f"<<<{final_answer}>>>")