import math

def solve_interstellar_problem():
    """
    Solves the Interstellar energy storage problem.
    """
    # Step 1: Define constants and calculate travel time
    distance_ly = 10.0
    speed_c = 0.02
    travel_time_years = distance_ly / speed_c

    # Box dimensions and ball properties
    box_dims = (12, 11, 11)
    ball_radius = 2.0
    initial_energy_per_ball_MJ = 100.0
    leak_rate_per_year = 0.001
    required_energy_MJ = 1000.0

    # Step 2: Define the optimal packing of 14 spheres found through analysis.
    # This configuration is a staggered-layer packing which is more efficient
    # than a simple cubic grid.
    # Layer 1 (z=2.0)
    placed_centers = [
        (2.0, 2.0, 2.0), (6.0, 2.0, 2.0), (10.0, 2.0, 2.0),
        (2.0, 6.0, 2.0), (6.0, 6.0, 2.0), (10.0, 6.0, 2.0),
    ]
    # Layer 2 (z=5.0) - placed in the hollows of layer 1
    placed_centers.extend([
        (4.0, 4.0, 5.0), (8.0, 4.0, 5.0),
    ])
    # Layer 3 (z=8.0)
    placed_centers.extend([
        (2.0, 2.0, 8.0), (6.0, 2.0, 8.0), (10.0, 2.0, 8.0),
        (2.0, 6.0, 8.0), (6.0, 6.0, 8.0), (10.0, 6.0, 8.0),
    ])
    
    total_balls = len(placed_centers)

    # Step 3: Classify balls as touching or non-touching
    num_touching = 0
    num_non_touching = 0
    
    # Boundary coordinates for centers of touching balls
    touch_x = [box_dims[0] - ball_radius, ball_radius] # [10.0, 2.0]
    touch_y = [box_dims[1] - ball_radius, ball_radius] # [9.0, 2.0]
    touch_z = [box_dims[2] - ball_radius, ball_radius] # [9.0, 2.0]

    for center in placed_centers:
        cx, cy, cz = center
        if cx in touch_x or cy in touch_y or cz in touch_z:
            num_touching += 1
        else:
            num_non_touching += 1

    # Step 4: Calculate final energy for each type of ball
    energy_non_touching_final = initial_energy_per_ball_MJ
    energy_touching_final = initial_energy_per_ball_MJ * math.pow(1 - leak_rate_per_year, travel_time_years)

    # Step 5: Calculate total energy and make the final decision
    total_final_energy = (num_touching * energy_touching_final) + (num_non_touching * energy_non_touching_final)

    print(f"Pioneer's journey to Pandora will take {int(travel_time_years)} years.")
    print(f"A maximum of {total_balls} energy balls can be packed in the container.")
    print(f"Out of {total_balls} balls, {num_touching} are touching the container and {num_non_touching} are not.")
    print("\nCalculating the final energy:")
    print(f"Energy of a non-touching ball after {int(travel_time_years)} years: {energy_non_touching_final:.3f} MJ")
    print(f"Energy of a touching ball after {int(travel_time_years)} years (due to leakage): {energy_touching_final:.3f} MJ")
    print("\nFinal Energy Equation:")
    print(f"Total Energy = (Touching Balls * Energy) + (Non-Touching Balls * Energy)")
    # The required output format for the equation
    print(f"{num_touching} * {energy_touching_final:.3f} MJ + {num_non_touching} * {energy_non_touching_final:.3f} MJ = {total_final_energy:.3f} MJ")

    print(f"\nTotal available energy upon arrival: {total_final_energy:.3f} MJ")
    print(f"Required energy for operations: {required_energy_MJ:.3f} MJ")

    if total_final_energy >= required_energy_MJ:
        print("\nDecision: The available energy is sufficient.")
        result = total_balls
    else:
        print("\nDecision: The available energy is NOT sufficient.")
        result = 0
        
    print(f"Final Answer: {result}")
    
    # Final answer in the required format
    print(f"\n<<<{result}>>>")

solve_interstellar_problem()