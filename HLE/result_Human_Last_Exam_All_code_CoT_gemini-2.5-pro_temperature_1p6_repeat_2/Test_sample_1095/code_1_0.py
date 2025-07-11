import sympy as sp

def solve_physics_problem():
    """
    This function explains the derivation for the condition on the radial wavevector k_r
    for a Bessel-Gauss beam to exhibit rotational propagation.
    """
    
    # Define symbols
    k_r, omega, beta, c, C1 = sp.symbols('k_r omega beta c C1')
    
    # Step 1: Explain the condition for rigid rotation without dispersion
    print("Step 1: To achieve rigid rotation of the wave packet, the group velocity of all constituent modes must be constant.")
    print("A sufficient condition for this is that the phase velocity (v_p) is constant for all modes.")
    print("v_p = omega / beta = constant (C1)\n")

    # Step 2: Use the dispersion relation for Bessel beams
    # beta = sqrt( (omega/c)^2 - k_r^2 )
    print("Step 2: The propagation constant (beta) for a Bessel beam is related to its frequency (omega) and radial wavevector (k_r) by:")
    print("beta = sqrt( (omega/c)^2 - k_r^2 )\n")

    # Step 3: Substitute and solve for k_r
    print("Step 3: Substitute beta into the constant phase velocity condition and solve for k_r.")
    # v_p = omega / sqrt((omega/c)**2 - k_r**2) = C1
    eq1 = sp.Eq(omega / sp.sqrt((omega/c)**2 - k_r**2), C1)
    print(f"Equation: {eq1}")
    
    # Square both sides
    eq2 = sp.Eq(eq1.lhs**2, eq1.rhs**2)
    # Solve for k_r**2
    k_r_squared = sp.solve(eq2, k_r**2)[0]
    
    print("Solving for k_r^2, we get:")
    print(f"k_r**2 = {k_r_squared}")
    print("This can be written as: k_r**2 = omega**2 * (1/c**2 - 1/C1**2)")
    print("Since (1/c**2 - 1/C1**2) is a constant, this means k_r**2 is proportional to omega**2.\n")
    
    # Step 4: Final Proportionality
    print("Step 4: Taking the square root, we find the proportionality for k_r.")
    print("k_r is proportional to omega (k_r ∝ ω).\n")

    # Step 5: Relate to topological charge l
    print("Step 5: In a 'light spring', frequency (omega) is linearly coupled to topological charge (l).")
    print("omega_l = omega_0 + l * delta_omega")
    print("This means omega is linearly proportional to l.\n")

    # Step 6: Conclusion
    print("Step 6: Combining these findings, since k_r ∝ omega and omega ∝ l, we conclude that:")
    print("k_r ∝ l")
    print("\nThe radial wavevector must be proportional to the topological charge of the beam.")

solve_physics_problem()