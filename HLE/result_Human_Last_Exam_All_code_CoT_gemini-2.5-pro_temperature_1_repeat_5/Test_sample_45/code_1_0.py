import sympy

def solve_telescope_problem():
    """
    This function derives the relationship between the focal length (f) of a
    liquid-mirror telescope and time (t), given a constant power source.
    It then determines the exponent 'n' in the proportionality f ∝ t^n.
    """

    # --- Step 1: Relate Focal Length (f) to Angular Velocity (ω) ---
    print("--- Step 1: Relate Focal Length (f) to Angular Velocity (ω) ---")
    print("The surface of a rotating liquid forms a paraboloid. Its height 'z' at a radial distance 'x' is given by:")
    print("z = (ω^2 * x^2) / (2 * g), where ω is the angular velocity and g is the acceleration due to gravity.")
    print("\nThe equation for a standard parabolic mirror with its vertex at the origin is x^2 = 4*f*z, where 'f' is the focal length.")
    print("Rearranging our liquid surface equation, we get: x^2 = (2*g / ω^2) * z.")
    print("By comparing the two equations for x^2, we can equate the coefficients of z:")
    print("4*f = 2*g / ω^2")
    print("Solving for f, we get: f = g / (2 * ω^2).")
    print("Since g is a constant, the focal length is inversely proportional to the square of the angular velocity.")
    print("This gives us our first key relationship: f ∝ (ω)^(-2)")

    # --- Step 2: Relate Angular Velocity (ω) to Time (t) ---
    print("\n--- Step 2: Relate Angular Velocity (ω) to Time (t) for a Constant Power Source ---")
    print("A constant power source 'P' is related to torque 'τ' and angular velocity 'ω' by: P = τ * ω.")
    print("Torque is also defined as τ = I * α, where 'I' is the moment of inertia and 'α' is the angular acceleration (dω/dt).")
    print("Substituting the expression for τ into the power equation gives: P = (I * dω/dt) * ω.")
    print("This is a separable differential equation. We rearrange it to solve for ω as a function of time t:")
    print("P * dt = I * ω * dω")
    print("We integrate both sides. The rotation starts from rest (ω=0 at t=0):")
    print("∫(P dt) from 0 to t  =  ∫(I * ω dω) from 0 to ω")
    print("The result of the integration is: P * t = (1/2) * I * ω^2.")
    print("Solving for ω^2, we find: ω^2 = (2 * P / I) * t.")
    print("Since P and I are constants, the square of the angular velocity is directly proportional to time.")
    print("This gives us our second key relationship: ω^2 ∝ t^1")

    # --- Step 3: Combine the Relationships to Find f(t) ---
    print("\n--- Step 3: Combine the Relationships to Find the dependence of f on t ---")
    print("We have two proportionalities:")
    print("1. f ∝ ω^(-2)")
    print("2. ω^2 ∝ t^1")
    print("We can rewrite the first proportionality as f ∝ (ω^2)^(-1).")
    print("Now, we substitute the second proportionality (ω^2 ∝ t^1) into the first:")
    print("f ∝ (t^1)^(-1)")
    print("This simplifies to: f ∝ t^(-1)")

    # --- Step 4: Determine the value of n ---
    print("\n--- Step 4: Determine the Final Value of n ---")
    print("The problem asks for the value of 'n' in the expression f ∝ t^n.")
    print("By comparing this with our derived relationship, we have the final equation:")
    
    # Using sympy to format the final output
    t = sympy.Symbol('t')
    f = sympy.Symbol('f')
    n = -1
    proportionality_eq = sympy.Eq(f, t**n)
    
    # Print the final equation with the value of n
    final_equation = f"f ∝ t^({n})"
    print(final_equation)
    print(f"\nBy comparison, the value of n is {n}.")

if __name__ == '__main__':
    solve_telescope_problem()
<<< -1 >>>