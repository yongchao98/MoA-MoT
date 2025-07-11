import sympy as sp

def solve_rod_sliding_angle():
    """
    This function derives and prints the expression for the angle at which a
    rotating rod begins to slide off a table edge.
    """
    # Define the symbolic variables used in the problem.
    # L: total length of the rod
    # l: distance from the pivot to the center of mass
    # mu: coefficient of static friction
    # theta: angle of the rod with the table
    L, l, mu = sp.symbols('L \u2113 \u03BC', positive=True, real=True)
    theta = sp.Symbol('\u03B8')

    # The derivation from the principles of rotational and translational dynamics
    # leads to the following relationship for tan(theta).
    #
    # The derivation steps involve:
    # 1. Moment of Inertia about pivot: I = M*(L**2/12 + l**2)
    # 2. Angular acceleration: alpha = M*g*l*cos(theta) / I
    # 3. Angular velocity squared: omega_sq = 2*g*l*sin(theta) / (L**2/12 + l**2)
    # 4. Normal Force: N = M*g*cos(theta) + M*l*omega_sq
    # 5. Friction Force: f = M*g*sin(theta) - M*l*alpha
    # 6. Sliding condition: f = mu * N
    #
    # Solving f = mu * N for tan(theta) yields the expression below.
    # We define an intermediate variable 'k' for clarity in the derivation.
    # k = l**2 / (l**2 + L**2 / 12)
    # tan_theta = (mu + k) / (1 - 2 * mu * k)
    #
    # Now, we construct this expression using sympy and simplify it.

    # Numerator of the expression for tan(theta) after substituting k and simplifying
    numerator = mu * L**2 + 12 * l**2 * (1 + mu)
    
    # Denominator of the expression
    denominator = L**2 + 12 * l**2 * (1 - 2 * mu)
    
    # Create the full expression for tan(theta)
    tan_theta_expression = numerator / denominator
    
    # Create the final equation object for printing
    final_equation = sp.Eq(sp.tan(theta), tan_theta_expression)

    # Print the final result in a readable format
    print("The angle \u03B8 at which the rod begins to slide is given by the equation:")
    sp.pprint(final_equation, use_unicode=True)

# Execute the function
solve_rod_sliding_angle()