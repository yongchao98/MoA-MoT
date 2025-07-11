import sympy

def calculate_max_distance():
    """
    Calculates and prints the symbolic expression for the maximum distance
    a spaceship can travel away from an asteroid, given a set of initial conditions
    and forces.
    """
    # Step 1: Define the symbolic variables for the problem.
    # G: Gravitational constant
    # m: mass of asteroid A
    # M: mass of spaceship B
    # l_0: initial distance between A and B
    # v_0: initial speed of spaceship B away from A
    # F: constant additional force applied to B in the direction of motion
    # l_max: the maximum distance to be calculated
    G, m, M, l_0, v_0, F = sympy.symbols('G m M l_0 v_0 F', real=True, positive=True)
    l_max = sympy.Symbol('l_max', real=True, positive=True)

    # Step 2: Use the Work-Energy Theorem to set up the governing equation.
    # The change in total energy is equal to the work done by the non-conservative force F.
    # E_total = Kinetic Energy + Gravitational Potential Energy
    # ΔE = W_F
    # (E_final - E_initial) = F * (l_max - l_0)
    # (0 - G*m*M/l_max) - (1/2*M*v_0**2 - G*m*M/l_0) = F*l_max - F*l_0
    #
    # Rearranging this gives a quadratic equation for l_max:
    # F*l_max**2 - (F*l_0 + G*m*M/l_0 - 1/2*M*v_0**2)*l_max + G*m*M = 0

    # Step 3: Define the coefficients of the quadratic equation a*x**2 + b*x + c = 0,
    # where x = l_max.
    a = F
    # For clarity, we define an intermediate constant 'C_term' which equals -b.
    C_term = F * l_0 + (G * m * M / l_0) - (sympy.Rational(1, 2) * M * v_0**2)
    b = -C_term
    c = G * m * M

    # Step 4: Solve for l_max using the quadratic formula.
    # The motion is from an initial point l_0 outwards to a turning point.
    # Therefore, l_max must be the larger of the two roots.
    solution = (-b + sympy.sqrt(b**2 - 4*a*c)) / (2*a)

    # Step 5: Print the components of the equation and the final solution.
    print("The maximum distance (l_max) can be found by solving a quadratic equation of the form:")
    print("a*(l_max)^2 + b*(l_max) + c = 0\n")
    
    print("The symbolic components of this equation are:")
    print("a = ", end="")
    sympy.pprint(a)
    print("\nb = ", end="")
    sympy.pprint(b)
    print("\nc = ", end="")
    sympy.pprint(c)
    
    print("\n-------------------------------------------------------------")
    print("The final expression for the maximum distance l_max is:")
    print("l_max = (-b + sqrt(b^2 - 4ac)) / (2a)")
    print("\nSubstituting the terms, we get:")
    sympy.pprint(solution)

calculate_max_distance()