def solve_focal_length_exponent():
    """
    This function derives and prints the relationship between the focal length 'f'
    of a rotating liquid mirror and time 't' when driven by a constant power source.
    It then determines the exponent 'n' in the relation f ∝ tⁿ.
    """

    print("--- Derivation of the focal length exponent 'n' ---")
    print("\nStep 1: Find the relationship between focal length (f) and angular speed (ω).")
    print("The surface of a liquid rotating with angular speed ω forms a paraboloid. The height 'h' of the surface at a radial distance 'x' from the axis of rotation is given by:")
    print("  h(x) = [ω² / (2g)] * x² + h₀")
    print("where 'g' is the acceleration due to gravity and h₀ is the height at the center.")
    print("\nThe equation for a standard parabolic mirror is y = a * x². Its focal length 'f' is given by f = 1 / (4a).")
    print("By comparing the two equations, we can identify a = ω² / (2g).")
    print("Substituting 'a' into the focal length formula:")
    print("  f = 1 / (4 * [ω² / (2g)]) = g / (2ω²)")
    print("This shows that the focal length is inversely proportional to the square of the angular speed: f ∝ ω⁻²")

    print("\n" + "-"*50 + "\n")

    print("Step 2: Find the relationship between angular speed (ω) and time (t).")
    print("The system is driven by a constant power source 'P' starting from rest.")
    print("Power is the rate of change of rotational kinetic energy (K): P = dK/dt.")
    print("The rotational kinetic energy is K = (1/2) * I * ω², where 'I' is the moment of inertia.")
    print("Since power 'P' is constant and the rotation starts from rest (K(0)=0), we can integrate with respect to time:")
    print("  ∫ P dt = ∫ dK  =>  P * t = K(t) - K(0) = K(t)")
    print("So, P * t = (1/2) * I * ω².")
    print("Solving for ω²:")
    print("  ω² = (2P / I) * t")
    print("This shows that ω² is directly proportional to time: ω² ∝ t")

    print("\n" + "-"*50 + "\n")

    print("Step 3: Combine the relationships to find f in terms of t.")
    print("From Step 1, we have: f ∝ ω⁻²")
    print("From Step 2, we have: ω² ∝ t")
    print("Substitute the expression for ω² into the proportionality for f:")
    print("  f ∝ (ω²)⁻¹ ∝ (t)⁻¹")
    print("So, the final relationship is: f ∝ t⁻¹")

    print("\n" + "-"*50 + "\n")
    
    # The exponent n in the relation f ∝ tⁿ is -1.
    n = -1

    print("Conclusion: Final Answer")
    print("We were asked to find the value of n in the relationship f ∝ tⁿ.")
    print(f"Based on the derivation, the final equation showing the proportionality is:")
    print(f"  f ∝ t^({n})")
    print(f"\nTherefore, the value of n is {n}.")

# Run the function to display the solution.
solve_focal_length_exponent()