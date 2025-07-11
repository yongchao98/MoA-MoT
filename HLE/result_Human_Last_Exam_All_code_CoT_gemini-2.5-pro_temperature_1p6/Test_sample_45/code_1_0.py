def derive_focal_length_exponent():
    """
    This script outlines the derivation to find the exponent 'n' in the relationship f ∝ t^n
    for a liquid-mirror telescope driven by a constant power source.
    The derivation is presented step-by-step.
    """

    print("Derivation of the focal length's time dependence for a liquid-mirror telescope.")
    print("==========================================================================")

    # --- Part 1: Relate focal length (f) to angular velocity (ω) ---
    print("\nStep 1: Express focal length (f) in terms of angular velocity (ω).")
    print("The surface of a rotating liquid forms a paraboloid given by the equation:")
    print("  z = (ω^2 * x^2) / (2 * g)")
    print("where 'z' is the height, 'x' is the radial distance, 'ω' is the angular velocity, and 'g' is the acceleration due to gravity.")
    print("\nThe equation for a parabolic mirror with focal length 'f' is:")
    print("  x^2 = 4 * f * z")
    print("\nBy comparing these two equations, we can find their relationship. Let's solve for z in the second equation and substitute into the first:")
    print("  z = x^2 / (4 * f)")
    print("  x^2 / (4 * f) = (ω^2 * x^2) / (2 * g)")
    print("Assuming x is not zero, we can simplify this to:")
    print("  1 / (4 * f) = ω^2 / (2 * g)")
    print("\nSolving for f, we find the relationship:")
    print("  f = g / (2 * ω^2)")
    print("This shows that the focal length 'f' is inversely proportional to the square of the angular velocity 'ω' (f ∝ ω^-2).")

    # --- Part 2: Relate angular velocity (ω) to time (t) ---
    print("\n--------------------------------------------------------------------------")
    print("\nStep 2: Express angular velocity (ω) in terms of time (t).")
    print("The system is driven by a constant power source, P.")
    print("According to the work-energy theorem, the work done on the system (W) equals its change in kinetic energy (ΔKE).")
    print("For a constant power source starting from rest at t=0, the work done is W = P * t.")
    print("The rotational kinetic energy of the fluid is KE = (1/2) * I * ω^2, where 'I' is the moment of inertia.")
    print("\nEquating work and kinetic energy:")
    print("  P * t = (1/2) * I * ω^2")
    print("\nSolving for ω^2, we get:")
    print("  ω^2 = (2 * P * t) / I")
    print("This shows that the square of the angular velocity, ω^2, is directly proportional to time 't' (ω^2 ∝ t).")

    # --- Part 3: Combine relationships and find the exponent n ---
    print("\n--------------------------------------------------------------------------")
    print("\nStep 3: Combine the results to find how f depends on t.")
    print("From Step 1, we have: f = g / (2 * ω^2)")
    print("From Step 2, we have: ω^2 = (2 * P * t) / I")
    print("\nSubstituting the expression for ω^2 into the equation for f:")
    print("  f = g / ( 2 * [ (2 * P * t) / I ] )")
    print("\nSimplifying the expression gives:")
    print("  f = (g * I) / (4 * P * t)")
    print("\nSince g, I, and P are constants, we can write the proportionality:")
    print("  f ∝ 1/t")

    # --- Part 4: Final Conclusion ---
    print("\n--------------------------------------------------------------------------")
    print("\nStep 4: Determine the value of n.")
    print("The problem states that the relationship is of the form: f ∝ t^n")
    print("Our derived relationship is: f ∝ 1/t, which is equivalent to f ∝ t^(-1).")

    n = -1
    print("\nBy comparing the exponents, we find the final equation with its number:")
    print(f"  f ∝ t^({n})")
    print("\nTherefore, the value of n is:")
    print(n)

# Execute the derivation function
derive_focal_length_exponent()