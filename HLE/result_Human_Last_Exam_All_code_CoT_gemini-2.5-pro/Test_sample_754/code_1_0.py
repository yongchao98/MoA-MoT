import math

def solve_interstellar_problem():
    """
    Calculates if the container box is sufficient for the Pioneer probe's mission to Pandora.
    """
    # Step 1: Define constants and calculate travel time
    distance_ly = 10.0  # light-years
    speed_c = 0.02      # proportion of speed of light
    travel_time_years = distance_ly / speed_c

    # Step 2: Define energy properties and calculate energy of a leaking ball
    initial_energy_per_ball_mj = 100.0
    leak_rate_per_year = 0.001
    
    # Energy of a ball touching the container wall after the journey
    energy_leaking_ball_mj = initial_energy_per_ball_mj * (1 - leak_rate_per_year)**travel_time_years
    
    # Energy of an inner ball remains the same
    energy_inner_ball_mj = initial_energy_per_ball_mj

    # Step 3: Define container and ball geometry and calculate max balls
    box_dims = {'x': 12.0, 'y': 11.0, 'z': 11.0}
    ball_radius = 2.0
    ball_diameter = 2.0 * ball_radius

    # Calculate how many balls fit in a simple grid packing
    num_balls_x = math.floor((box_dims['x'] - 2 * ball_radius) / ball_diameter) + 1
    num_balls_y = math.floor((box_dims['y'] - 2 * ball_radius) / ball_diameter) + 1
    num_balls_z = math.floor((box_dims['z'] - 2 * ball_radius) / ball_diameter) + 1
    
    total_balls = num_balls_x * num_balls_y * num_balls_z

    # Step 4: Classify balls as "touching" or "inner"
    # An inner ball is not on the edge of the packed grid in any dimension.
    # The grid of balls has dimensions num_balls_x, num_balls_y, num_balls_z.
    # An inner core exists only if all dimensions are > 2.
    if num_balls_x > 2 and num_balls_y > 2 and num_balls_z > 2:
        num_inner_balls = (num_balls_x - 2) * (num_balls_y - 2) * (num_balls_z - 2)
    else:
        # A more general approach to find inner balls:
        # We place ball centers starting at 'radius' with 'diameter' spacing.
        # An inner ball's center must not be at the first or last position on any axis.
        # Center x-coords: 2, 6, 10. Inner is {6}. (1 inner position)
        # Center y-coords: 2, 6. No inner positions. (0 inner positions)
        # Wait, the condition is not touching the *grid edge*, but the *box wall*.
        # A ball is touching if center_coord = radius OR center_coord = dim - radius
        # Let's find the single inner ball at the center of the box if it exists.
        x_center_pos = ball_radius + (num_balls_x // 2) * ball_diameter
        y_center_pos = ball_radius + (num_balls_y // 2) * ball_diameter
        z_center_pos = ball_radius + (num_balls_z // 2) * ball_diameter

        # A ball is inner if its center is not on a boundary plane
        is_inner = True
        if x_center_pos == ball_radius or x_center_pos == box_dims['x'] - ball_radius: is_inner = False
        if y_center_pos == ball_radius or y_center_pos == box_dims['y'] - ball_radius: is_inner = False
        if z_center_pos == ball_radius or z_center_pos == box_dims['z'] - ball_radius: is_inner = False
        
        # In our case, the packing grid is 3x2x2. Only the ball at (6,6,6) is truly inner.
        # Any ball with center coord y=2 or z=2 will touch the wall.
        # Any ball with center coord x=2 or x=10 will touch the wall.
        # The ball at (6,6,6) is the only one with no touching coordinates.
        num_inner_balls = 1
        
    num_touching_balls = total_balls - num_inner_balls

    # Step 5: Calculate total energy
    required_energy_mj = 1000.0
    total_energy_mj = (num_inner_balls * energy_inner_ball_mj) + \
                      (num_touching_balls * energy_leaking_ball_mj)

    # Step 6: Make decision and print output
    print("Pioneer Mission Energy Analysis:")
    print(f"Maximal balls that can be packed: {total_balls}")
    print(f"Number of inner balls: {num_inner_balls}")
    print(f"Number of touching balls: {num_touching_balls}")
    print(f"Energy of one inner ball on arrival: {energy_inner_ball_mj:.4f} MJ")
    print(f"Energy of one touching ball on arrival: {energy_leaking_ball_mj:.4f} MJ")
    print("\nFinal Energy Calculation:")
    # The problem asks to print the numbers in the final equation
    print(f"({num_inner_balls} * {energy_inner_ball_mj}) + ({num_touching_balls} * {energy_leaking_ball_mj}) = {total_energy_mj} MJ")

    if total_energy_mj >= required_energy_mj:
        print(f"\nTotal energy {total_energy_mj:.2f} MJ is SUFFICIENT (>= {required_energy_mj} MJ).")
        final_answer = total_balls
    else:
        print(f"\nTotal energy {total_energy_mj:.2f} MJ is NOT SUFFICIENT (< {required_energy_mj} MJ).")
        final_answer = 0
    
    # Final answer in the required format
    print(f"\n<<<{final_answer}>>>")

solve_interstellar_problem()