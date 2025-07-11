import math
import numpy as np

def solve_interstellar_problem():
    """
    Solves the Interstellar energy storage problem by calculating travel time,
    energy decay, and running a greedy packing algorithm to find the maximum
    achievable energy.
    """
    # 1. Define constants and calculate travel-related values
    TRAVEL_DISTANCE_LY = 10.0
    SPEED_FRAC_C = 0.02
    BOX_DIMS = (12, 11, 11)
    BALL_RADIUS = 2.0
    BALL_DIAMETER = 4.0
    COORD_STEP = 0.5
    INITIAL_ENERGY = 100.0  # MJ
    LEAK_RATE_PER_YEAR = 0.001
    REQUIRED_ENERGY = 1000.0  # MJ

    # Calculate time and energy of balls upon arrival
    travel_time_years = TRAVEL_DISTANCE_LY / SPEED_FRAC_C
    total_leak_fraction = LEAK_RATE_PER_YEAR * travel_time_years
    energy_leaking_ball = INITIAL_ENERGY * (1 - total_leak_fraction)
    energy_inner_ball = INITIAL_ENERGY

    # 2. Define the search space for ball centers
    min_cx = BALL_RADIUS
    max_cx = BOX_DIMS[0] - BALL_RADIUS
    min_cy = BALL_RADIUS
    max_cy = BOX_DIMS[1] - BALL_RADIUS
    min_cz = BALL_RADIUS
    max_cz = BOX_DIMS[2] - BALL_RADIUS

    # 3. Generate and sort all possible center coordinates for the greedy algorithm
    possible_centers = []
    x_coords = np.arange(min_cx, max_cx + COORD_STEP, COORD_STEP)
    y_coords = np.arange(min_cy, max_cy + COORD_STEP, COORD_STEP)
    z_coords = np.arange(min_cz, max_cz + COORD_STEP, COORD_STEP)

    for x in x_coords:
        for y in y_coords:
            for z in z_coords:
                possible_centers.append((x, y, z))

    # Sort candidates by distance from the geometric center to prioritize inner placements
    box_center = ( (min_cx + max_cx) / 2.0, (min_cy + max_cy) / 2.0, (min_cz + max_cz) / 2.0 )
    
    def dist_sq_from_center(p):
        return (p[0] - box_center[0])**2 + (p[1] - box_center[1])**2 + (p[2] - box_center[2])**2
        
    possible_centers.sort(key=dist_sq_from_center)

    # 4. Run the greedy packing algorithm
    placed_balls_centers = []
    min_dist_sq = BALL_DIAMETER**2

    for candidate_center in possible_centers:
        can_place = True
        for placed_center in placed_balls_centers:
            dist_sq = (candidate_center[0] - placed_center[0])**2 + \
                      (candidate_center[1] - placed_center[1])**2 + \
                      (candidate_center[2] - placed_center[2])**2
            if dist_sq < min_dist_sq:
                can_place = False
                break
        if can_place:
            placed_balls_centers.append(candidate_center)

    # 5. Classify balls and calculate total energy
    num_inner_balls = 0
    num_touching_balls = 0
    for x, y, z in placed_balls_centers:
        # A ball is touching if its center is at the boundary of the allowed space for centers
        is_touching = (
            math.isclose(x, min_cx) or math.isclose(x, max_cx) or
            math.isclose(y, min_cy) or math.isclose(y, max_cy) or
            math.isclose(z, min_cz) or math.isclose(z, max_cz)
        )
        if is_touching:
            num_touching_balls += 1
        else:
            num_inner_balls += 1
    
    total_energy_achieved = (num_inner_balls * energy_inner_ball) + (num_touching_balls * energy_leaking_ball)
    
    # Print the detailed breakdown of the calculation
    print("Final Energy Calculation:")
    print(f"Travel time = {TRAVEL_DISTANCE_LY:.0f} light-years / {SPEED_FRAC_C}c = {travel_time_years:.0f} years")
    print(f"Energy of a leaking ball = {INITIAL_ENERGY:.1f} MJ * (1 - {LEAK_RATE_PER_YEAR} * {travel_time_years:.0f}) = {energy_leaking_ball:.1f} MJ")
    print(f"Energy of an inner ball = {energy_inner_ball:.1f} MJ")
    print(f"Greedy packing resulted in {len(placed_balls_centers)} balls:")
    print(f"  - Inner balls (no leak): {num_inner_balls}")
    print(f"  - Touching balls (leaking): {num_touching_balls}")
    print(f"Total energy = {num_inner_balls} * {energy_inner_ball:.1f} + {num_touching_balls} * {energy_leaking_ball:.1f} = {total_energy_achieved:.1f} MJ")
    print(f"Required energy for operations = {REQUIRED_ENERGY:.1f} MJ")

    # 6. Determine the final answer
    if total_energy_achieved >= REQUIRED_ENERGY:
        print("Conclusion: The stored energy is sufficient for operations.")
        return len(placed_balls_centers)
    else:
        print("Conclusion: The stored energy is NOT sufficient for operations.")
        return 0

# Run the solver and print the final answer in the required format.
final_answer = solve_interstellar_problem()
print(f"\n<<<{final_answer}>>>")
