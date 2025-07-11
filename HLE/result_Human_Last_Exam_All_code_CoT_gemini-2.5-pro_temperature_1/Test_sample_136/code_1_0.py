import math

def solve_torus_walk_problem():
    """
    Calculates the limit of the conditional probability for a random walk on a 2D torus.

    The problem asks for the limit as n -> infinity of:
    P[x_0 was not visited before t_n | 0 was not visited before t_n]
    where t_n = n^2 * (ln n)^2.

    The solution uses concepts from potential theory of random walks, namely
    capacity and its relation to electrical resistance.
    """

    # The conditional probability P_n can be expressed asymptotically as:
    # P_n ~ exp(-t_n * (Cap({0, x_0}) - Cap({0})))

    # The difference in capacity is given by:
    # delta_Cap = pi(x_0) * P_x0(T_0 > T_x0^+)
    # where pi(x_0) = 1/n^2 for the torus.
    # So, delta_Cap = (1/n^2) * (1 - P_x0(T_0 < T_x0^+))
    # The exponent is -t_n * delta_Cap = -(ln n)^2 * (1 - P_x0(T_0 < T_x0^+))

    # The hitting probability is related to electrical resistance R(0, x_0):
    # P_x0(T_0 < T_x0^+) = 1 / (deg(x_0) * R(0, x_0))
    # The degree of any vertex in the 2D torus is 4.

    # The resistance R(0, x_0) for x_0=(1,1) on a large torus is
    # approximated by the resistance on an infinite 2D grid, which is 2/pi.
    
    deg_x0 = 4
    R_0_x0_val = 2 / math.pi
    
    # Calculate the hitting probability
    p_hit_before_return = 1 / (deg_x0 * R_0_x0_val)
    
    # The final expression for the exponent is -(ln n)^2 * C
    C = 1 - p_hit_before_return

    print("The conditional probability P_n has the asymptotic form: exp(- (ln n)^2 * C)")
    print("The constant C is derived as follows:")
    print(f"C = 1 - P_x0(T_0 < T_x0^+) = 1 - 1 / (deg(x_0) * R(0, x_0))")
    print(f"Using deg(x_0) = {deg_x0} and R(0, x_0) ≈ 2/π:")
    print(f"C = 1 - 1 / ({deg_x0} * (2/{math.pi:.4f}))")
    print(f"C = 1 - π/8")
    print(f"Numerically, C ≈ {C:.4f}")

    # As n -> infinity, (ln n)^2 -> infinity.
    # Since C is a positive constant, the exponent -(ln n)^2 * C -> -infinity.
    # The limit of the exponential of a value going to -infinity is 0.
    final_limit = 0
    
    print("\nSince C is positive, the exponent -(ln n)^2 * C approaches -infinity.")
    print("Therefore, the limit of the probability is 0.")
    print("\nThe final equation for the probability P_n is:")
    print(f"P_n ≈ exp( - (ln n)^2 * (1 - π/8) )")
    print("\nThe numbers in this equation are 1 and 8, and the constant π.")
    print(f"\nThe final limit is: {final_limit}")

solve_torus_walk_problem()
<<<0>>>