import math

def solve_interstellar_problem():
    """
    Calculates if the Pioneer probe's energy container is sufficient for its mission to Pandora.
    """
    # --- Mission and Environment Parameters ---
    distance_ly = 10.0
    speed_frac_c = 0.02
    required_energy_mj = 1000.0

    # --- Energy Ball and Container Parameters ---
    ball_radius_cm = 2.0
    ball_diameter_cm = ball_radius_cm * 2
    initial_ball_energy_mj = 100.0
    leak_rate_per_year = 0.001  # 0.1%
    box_dims_cm = (12.0, 11.0, 11.0)  # (Length, Width, Height)

    # Step 1: Calculate travel time in years
    travel_time_years = distance_ly / speed_frac_c
    print(f"The probe will travel for {travel_time_years:.0f} years.")

    # Step 2: Calculate maximum number of balls that can be packed
    box_l, box_w, box_h = box_dims_cm
    
    # The space available for the centers of the balls
    space_l = box_l - 2 * ball_radius_cm
    space_w = box_w - 2 * ball_radius_cm
    space_h = box_h - 2 * ball_radius_cm
    
    # Max balls that fit along each dimension in a grid
    n_l = math.floor(space_l / ball_diameter_cm) + 1
    n_w = math.floor(space_w / ball_diameter_cm) + 1
    n_h = math.floor(space_h / ball_diameter_cm) + 1
    n_total = n_l * n_w * n_h
    
    print(f"The container can hold a grid of {n_l}x{n_w}x{n_h} = {n_total} balls.")

    # Step 3: Differentiate between surface (leaking) and interior (non-leaking) balls
    
    # A ball is interior along an axis if the grid of balls doesn't fill the entire available center space.
    # Slack space = (available center space) - (grid extent)
    slack_l = space_l - (n_l - 1) * ball_diameter_cm
    slack_w = space_w - (n_w - 1) * ball_diameter_cm
    slack_h = space_h - (n_h - 1) * ball_diameter_cm

    # If slack space is 0, only the balls not on the ends of the grid are interior along that axis.
    # If slack space is > 0, the entire grid can be centered, making all balls interior along that axis.
    n_interior_layers_l = n_l if slack_l > 0.0001 else max(0, n_l - 2) # Use a small tolerance for float comparison
    n_interior_layers_w = n_w if slack_w > 0.0001 else max(0, n_w - 2)
    n_interior_layers_h = n_h if slack_h > 0.0001 else max(0, n_h - 2)
    
    n_interior = n_interior_layers_l * n_interior_layers_w * n_interior_layers_h
    n_surface = n_total - n_interior
    
    print(f"Number of interior (non-leaking) balls: {n_interior}")
    print(f"Number of surface (leaking) balls: {n_surface}")

    # Step 4: Calculate the total energy upon arrival
    energy_leaking_ball = initial_ball_energy_mj * math.pow(1 - leak_rate_per_year, travel_time_years)
    energy_interior_ball = initial_ball_energy_mj

    total_energy = (n_surface * energy_leaking_ball) + (n_interior * energy_interior_ball)
    
    print(f"\nCalculating the final total energy:")
    print(f"Energy of one interior ball at arrival: {energy_interior_ball:.2f} MJ")
    print(f"Energy of one surface ball at arrival: {energy_leaking_ball:.2f} MJ")
    
    # Print the equation with numbers
    term1 = n_surface * energy_leaking_ball
    term2 = n_interior * energy_interior_ball
    print(f"Total Energy = ({n_surface} * {energy_leaking_ball:.2f}) + ({n_interior} * {energy_interior_ball:.2f})")
    print(f"Total Energy = {term1:.2f} MJ + {term2:.2f} MJ = {total_energy:.2f} MJ")

    # Step 5: Final Verdict
    print(f"\nRequired energy for operations is {required_energy_mj:.2f} MJ.")
    
    if total_energy >= required_energy_mj:
        print(f"The total energy of {total_energy:.2f} MJ is sufficient.")
        final_answer = n_total
    else:
        print(f"The total energy of {total_energy:.2f} MJ is not sufficient.")
        final_answer = 0
        
    print(f"\nFinal Answer: {final_answer}")
    print(f"<<<{final_answer}>>>")

if __name__ == '__main__':
    solve_interstellar_problem()