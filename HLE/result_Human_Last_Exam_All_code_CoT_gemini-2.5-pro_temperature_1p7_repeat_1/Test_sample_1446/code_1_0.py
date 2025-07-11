def solve_critical_exponent_nu():
    """
    Calculates the critical exponent nu (ν) within the mean-field approximation
    of the Ginzburg-Landau (φ⁴) theory.

    This value is exact for spatial dimensions d >= 4.
    """

    print("Solving for the critical exponent ν in a G₄ (φ⁴) framework using Mean-Field Theory.")
    print("-" * 75)

    # Step 1: State the Ginzburg-Landau free energy density in momentum space (quadratic part).
    # F(q) ≈ r * φ(q)² + q² * φ(q)² = (r + q²) * φ(q)²
    # The propagator or correlation function G(q) is the inverse of the quadratic part.
    # G(q) ∝ 1 / (r + q²)
    print("Step 1: In mean-field theory, the correlation function G(q) in momentum space 'q' is proportional to 1 / (r + q²), where 'r' is the coefficient of the φ² term and is proportional to the reduced temperature (T - T_c).")

    # Step 2: Define the correlation length (ξ).
    # The standard form for the correlation function is G(q) ∝ 1 / (ξ⁻² + q²).
    # By comparing the two forms, we can relate the correlation length ξ to the parameter r.
    print("\nStep 2: The correlation length ξ is defined such that G(q) ∝ 1 / (ξ⁻² + q²).")
    print("Comparing the two forms gives the relationship: ξ⁻² = r.")

    # Step 3: Relate ξ to the reduced temperature t = |T - T_c|.
    # From mean-field theory, we know r is directly proportional to the reduced temperature.
    # r ∝ |T - T_c|¹
    # Therefore, ξ⁻² ∝ |T - T_c|¹
    power_of_t = 1
    print(f"\nStep 3: From the theory, we know that r ∝ |T - T_c|^{power_of_t}.")
    print(f"Therefore, we have the relation: ξ⁻² ∝ |T - T_c|^{power_of_t}.")
    
    # Step 4: Solve for ξ.
    # (ξ²)⁻¹ ∝ |T - T_c|¹  =>  ξ² ∝ |T - T_c|⁻¹  => ξ ∝ (|T - T_c|⁻¹)¹/²
    # ξ ∝ |T - T_c|⁻¹/²
    exponent_of_xi = -power_of_t / 2.0
    print(f"\nStep 4: Solving for ξ gives: ξ ∝ |T - T_c|^({exponent_of_xi}).")

    # Step 5: Compare with the definition of the critical exponent ν.
    # The definition of ν is: ξ ∝ |T - T_c|⁻ν
    # Comparing exponents, we get: -ν = -1/2
    # So, ν = 1/2
    print("\nStep 5: The critical exponent ν is defined by the scaling relation: ξ ∝ |T - T_c|⁻ν.")
    print("By comparing the derived exponent with the definition, we find -ν = -1/2.")

    # Final Result
    numerator = 1
    denominator = 2
    nu_value = numerator / denominator

    print("\n" + "-" * 75)
    print("Final Calculation:")
    print(f"The equation for the exponent is: ν = {numerator} / {denominator}")
    print(f"The value of the critical exponent ν is: {nu_value}")
    print("-" * 75)


if __name__ == "__main__":
    solve_critical_exponent_nu()
