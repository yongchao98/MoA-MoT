import math

def solve_parliament_design():
    """
    Calculates the maximum integer value for the parameter K in the parliament design.
    """
    # Step 1 & 2: Define radial positions
    # Number of rows is 13 (791 members / 61 sections)
    # r_n = 3 + (n-1) * 1.5
    r_1 = 3.0  # Speaker's row
    r_13 = 3.0 + (13 - 1) * 1.5  # Observer's row

    # Step 5 & 6: Find the most restrictive constraint and solve for K
    # The visibility constraint must hold for all intermediate rows (n=2 to 12).
    # The inequality is K < 2 * (r_13 - r_1) * (r_n - r_1).
    # The most restrictive (lowest) upper bound for K occurs when (r_n - r_1) is minimized.
    # This happens for the obstacle closest to the speaker, which is the person in row 2 (n=2).
    
    r_2 = 3.0 + (2 - 1) * 1.5 # Obstacle's row
    
    # Calculate the components of the inequality
    r_observer_minus_r_speaker = r_13 - r_1
    r_obstacle_minus_r_speaker = r_2 - r_1
    
    # Calculate the upper bound for K
    K_upper_bound = 2 * r_observer_minus_r_speaker * r_obstacle_minus_r_speaker
    
    # The maximum integer value K can take is the floor of this upper bound.
    max_integer_K = math.floor(K_upper_bound)

    # Output the final equation with numbers
    print("The visibility constraint is given by the inequality: K < 2 * (r_observer - r_speaker) * (r_obstacle - r_speaker)")
    print(f"Substituting the values for the most restrictive case (obstacle in row 2):")
    print(f"K < 2 * ({r_13} - {r_1}) * ({r_2} - {r_1})")
    print(f"K < 2 * {r_observer_minus_r_speaker} * {r_obstacle_minus_r_speaker}")
    print(f"K < {K_upper_bound}")
    print(f"\nThe maximum integer value K can take is {max_integer_K}.")

solve_parliament_design()
<<<53>>>