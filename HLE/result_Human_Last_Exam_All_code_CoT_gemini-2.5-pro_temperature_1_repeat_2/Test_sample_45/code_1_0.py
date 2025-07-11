import sympy as sp

def solve_telescope_focal_length_exponent():
    """
    This function symbolically derives the relationship between focal length f and time t
    for a liquid-mirror telescope spun up by a constant power source.
    It then computes and prints the exponent n in the relation f ∝ t^n.
    """

    print("Deriving the exponent 'n' for the focal length relationship f ∝ t^n.\n")

    # --- Define symbolic variables ---
    # t: time, P: constant power, I: moment of inertia
    # g: acceleration due to gravity, f: focal length
    # omega_sq: square of the angular velocity (ω^2)
    t, P, I, g = sp.symbols('t P I g', positive=True)
    f_symbol = sp.Symbol('f')
    omega_sq = sp.Symbol('omega_sq', positive=True)

    # --- Step 1: Relate power, time, and angular velocity ---
    # The kinetic energy K = (1/2) * I * ω^2.
    # For constant power P, K = P * t (integrating dK/dt = P from rest).
    # This gives the equation: (1/2) * I * ω^2 = P * t.
    print("Step 1: Relate kinetic energy (K = P*t) to rotational energy (K = 1/2 * I * ω^2).")
    kinetic_energy_eq = sp.Eq(sp.Rational(1, 2) * I * omega_sq, P * t)
    print(f"Equation: {kinetic_energy_eq}")
    
    # Solve for omega_sq (ω^2) in terms of time t
    omega_sq_solution = sp.solve(kinetic_energy_eq, omega_sq)[0]
    print(f"Solving for ω^2 gives: ω^2 = {omega_sq_solution}\n")

    # --- Step 2: Relate focal length and angular velocity ---
    # The surface of a rotating liquid is a paraboloid z = (ω^2 / (2*g)) * x^2.
    # The focal length f of such a parabola is f = g / (2 * ω^2).
    print("Step 2: Relate focal length (f) to the square of the angular velocity (ω^2).")
    focal_length_relation = g / (2 * omega_sq)
    print(f"Equation: f = {focal_length_relation}\n")

    # --- Step 3: Find focal length as a function of time ---
    # Substitute the expression for ω^2 from Step 1 into the focal length equation from Step 2.
    print("Step 3: Substitute ω^2 from Step 1 into the focal length equation from Step 2.")
    f_vs_t = focal_length_relation.subs(omega_sq, omega_sq_solution)
    final_equation = sp.Eq(f_symbol, f_vs_t)
    print("The final equation for focal length f as a function of time t is:")
    sp.pretty_print(final_equation)
    print("")

    # --- Step 4: Determine the exponent n ---
    # The relationship is of the form f ∝ t^n.
    # We compute n formally using the property that n = (t/f) * (df/dt).
    print("Step 4: Determine the exponent 'n' from the relationship f ∝ t^n.")
    n = sp.simplify((t / f_vs_t) * sp.diff(f_vs_t, t))
    
    print("The final relationship is f ∝ t^n. The number in the final equation is the exponent n.")
    print("The computed value for n is:")
    print(int(n))

# Execute the function to find and print the exponent n.
solve_telescope_focal_length_exponent()