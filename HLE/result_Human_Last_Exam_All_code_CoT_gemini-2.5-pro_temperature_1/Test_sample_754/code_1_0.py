import math

def solve_interstellar_problem():
    """
    Calculates if the energy container is sufficient for the Pioneer probe's mission.
    """
    # Problem parameters
    distance_ly = 10.0
    speed_fraction_c = 0.02
    required_energy_mj = 1000.0
    initial_energy_per_ball_mj = 100.0
    leak_rate_per_year = 0.001
    ball_radius_cm = 2.0
    box_dims_cm = (12.0, 11.0, 11.0)
    
    # --- Step 1: Calculate travel time ---
    travel_time_years = distance_ly / speed_fraction_c

    # --- Step 2: Calculate maximum number of balls and their positions ---
    def get_center_coords(dimension, radius):
        """Calculates possible center coordinates along one dimension."""
        coords = []
        pos = radius
        while pos <= dimension - radius:
            coords.append(pos)
            pos += 2 * radius
        return coords

    x_coords = get_center_coords(box_dims_cm[0], ball_radius_cm)
    y_coords = get_center_coords(box_dims_cm[1], ball_radius_cm)
    z_coords = get_center_coords(box_dims_cm[2], ball_radius_cm)
    
    max_balls = len(x_coords) * len(y_coords) * len(z_coords)

    # --- Step 3: Differentiate between touching and non-touching balls ---
    num_non_touching = 0
    for x in x_coords:
        for y in y_coords:
            for z in z_coords:
                is_touching_x = (x - ball_radius_cm <= 0) or (x + ball_radius_cm >= box_dims_cm[0])
                is_touching_y = (y - ball_radius_cm <= 0) or (y + ball_radius_cm >= box_dims_cm[1])
                is_touching_z = (z - ball_radius_cm <= 0) or (z + ball_radius_cm >= box_dims_cm[2])
                
                if not is_touching_x and not is_touching_y and not is_touching_z:
                    num_non_touching += 1
    
    num_touching = max_balls - num_non_touching

    # --- Step 4: Calculate the energy remaining after the journey ---
    energy_non_touching_ball = initial_energy_per_ball_mj
    
    # Calculate energy for a single touching ball after decay
    remaining_factor = (1 - leak_rate_per_year) ** travel_time_years
    energy_touching_ball = initial_energy_per_ball_mj * remaining_factor

    # Calculate total energy
    total_energy_mj = (num_non_touching * energy_non_touching_ball) + (num_touching * energy_touching_ball)

    # --- Step 5: Final conclusion and output ---
    print(f"Analysis for Pioneer's journey to Pandora:")
    print("-" * 40)
    print(f"1. Travel Time: {travel_time_years:.0f} years")
    print(f"2. Max Balls: {max_balls} ({len(x_coords)}x{len(y_coords)}x{len(z_coords)})")
    print(f"   - Non-touching balls: {num_non_touching}")
    print(f"   - Touching balls: {num_touching}")
    print("-" * 40)
    
    print("Final Energy Calculation:")
    print(f"Energy per non-touching ball = {energy_non_touching_ball:.1f} MJ")
    print(f"Energy per touching ball = {initial_energy_per_ball_mj:.1f} * (1 - {leak_rate_per_year})**{travel_time_years:.0f} = {energy_touching_ball:.3f} MJ")
    
    print("\nTotal Energy on Arrival Equation:")
    # The user requested to see the final equation with all numbers
    print(f"E_total = ({num_non_touching} balls * {energy_non_touching_ball:.1f} MJ) + ({num_touching} balls * {energy_touching_ball:.3f} MJ)")
    print(f"E_total = {num_non_touching * energy_non_touching_ball:.3f} MJ + {num_touching * energy_touching_ball:.3f} MJ")
    print(f"E_total = {total_energy_mj:.3f} MJ")
    
    print("\nConclusion:")
    print(f"The total energy on arrival ({total_energy_mj:.3f} MJ) is less than the required {required_energy_mj:.1f} MJ.")
    
    if total_energy_mj >= required_energy_mj:
        final_answer = max_balls
    else:
        final_answer = 0
    
    # The final answer in the required format
    print(f"\nIs the container sufficient? -> No.")
    print(f"<<<{final_answer}>>>")

solve_interstellar_problem()