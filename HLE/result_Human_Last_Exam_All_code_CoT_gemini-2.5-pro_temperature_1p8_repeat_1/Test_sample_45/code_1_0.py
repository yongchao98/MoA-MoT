def solve_telescope_problem():
    """
    This program outlines the derivation to find the exponent 'n' in the
    focal length relationship f ∝ t^n for a liquid-mirror telescope.
    """
    
    print("Step 1: Find the relationship between focal length (f) and angular velocity (ω).")
    print("The surface of the rotating liquid is a parabola defined by: z = (ω² * r²) / (2g)")
    print("The standard equation for a parabolic mirror with focal length f is: r² = 4 * f * z")
    print("By comparing these two equations, we set their terms for r² equal:")
    print("4 * f * z = (2g / ω²) * z")
    print("Solving for f, we get: f = 2g / (4 * ω²) = g / (2 * ω²)")
    print("This shows that f is proportional to ω⁻².\n")

    print("Step 2: Find the relationship between angular velocity (ω) and time (t).")
    print("A constant power (P) is applied. Power is the rate of change of kinetic energy (K).")
    print("Rotational kinetic energy K = (1/2) * I * ω², where I is the moment of inertia.")
    print("P = dK/dt = d/dt( (1/2) * I * ω² ) = I * ω * (dω/dt)")
    print("This gives the differential equation: I * ω * dω = P * dt")
    print("Integrating both sides (from t=0, ω=0 to t, ω): ∫(I * ω')dω' = ∫(P)dt'")
    print("The result of integration is: (1/2) * I * ω² = P * t")
    print("Solving for ω², we get: ω² = (2 * P * t) / I")
    print("This shows that ω² is proportional to t.\n")

    print("Step 3: Combine the results to find how f depends on t.")
    print("From Step 1, we have: f ∝ ω⁻²")
    print("From Step 2, we have: ω² ∝ t")
    print("Substituting the relation for ω² into the relation for f:")
    print("f ∝ (t)⁻¹\n")

    print("Conclusion:")
    print("The relationship between focal length and time is f ∝ tⁿ.")
    n = -1
    print(f"By comparing this with our derived relation f ∝ t⁻¹, we find the value of n.")
    
    # As requested, output the final equation showing all numbers.
    print(f"\nThe final equation is f ∝ t^({n}).")
    print(f"Therefore, n = {n}.")

solve_telescope_problem()