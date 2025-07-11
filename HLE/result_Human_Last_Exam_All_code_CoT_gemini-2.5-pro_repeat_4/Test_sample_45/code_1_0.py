import sympy

def solve_telescope_problem():
    """
    This function outlines the derivation for the focal length of a liquid-mirror
    telescope driven by a constant power source and prints the final result.
    """

    # Define symbols for clarity in the explanation
    f, g, omega, P, I, t = sympy.symbols('f g ω P I t')
    k1, k2 = sympy.symbols('k1 k2')

    # Step 1: Relationship between focal length (f) and angular velocity (ω)
    # The surface of a rotating liquid is a paraboloid: z = ω^2 * r^2 / (2*g)
    # The equation of a parabolic mirror is r^2 = 4*f*z.
    # Equating them gives f = g / (2*ω^2).
    print("Step 1: The focal length (f) is inversely proportional to the square of the angular velocity (ω).")
    print("f ∝ 1/ω^2")
    # This can be written as: f = k1 * ω^(-2), where k1 is a constant (g/2).
    print(f"Equation 1: f = {k1} * ω**(-2)\n")

    # Step 2: Relationship between angular velocity (ω) and time (t)
    # For a constant power source P, we have P = Torque * ω.
    # Torque = I * dω/dt, where I is the moment of inertia (assumed constant).
    # This leads to the differential equation: P = I * ω * dω/dt.
    # Integrating P*dt = I*ω*dω gives P*t = (1/2)*I*ω^2.
    # Solving for ω^2 gives: ω^2 = (2*P/I)*t.
    print("Step 2: For a constant power source, the square of the angular velocity (ω^2) is proportional to time (t).")
    print("ω^2 ∝ t")
    # This can be written as: ω^2 = k2 * t, where k2 is a constant (2*P/I).
    print(f"Equation 2: ω^2 = {k2} * t\n")

    # Step 3: Combine the relationships to find f(t)
    # Substitute Equation 2 into Equation 1.
    # f = k1 / (k2 * t) = (k1/k2) * t^(-1)
    print("Step 3: By substituting ω^2 from Step 2 into the equation from Step 1, we find the relationship between f and t.")
    print("f ∝ 1/t")
    print("f ∝ t^(-1)\n")

    # Step 4: Determine the exponent n
    # The problem asks for the exponent n in the relationship f ∝ t^n.
    # From our derivation, the exponent is -1.
    n = -1
    print("The final relationship is of the form: f ∝ t^n")
    print(f"Based on the derivation, the final equation is f ∝ t^({n}).")
    print(f"The value of the exponent n is: {n}")

solve_telescope_problem()