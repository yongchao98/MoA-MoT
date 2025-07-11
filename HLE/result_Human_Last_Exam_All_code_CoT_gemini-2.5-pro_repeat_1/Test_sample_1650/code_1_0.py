def solve_overlap_integral():
    """
    This function prints the final analytical expression for the overlap
    integral of two 2s orbitals in the H₂⁺ ion.
    """
    # Define the coefficients of the polynomial in the final expression
    c0 = 1
    c1 = 1/2
    c2 = 1/12
    c3 = 0  # There is no ρ³ term in the final result
    c4 = 1/240

    # Print an explanation of the formula
    print("The analytical expression for the overlap integral S between two 2s orbitals is:")
    print("S(ρ) = exp(-ρ/2) * (c₀ + c₁*ρ + c₂*ρ² + c₃*ρ³ + c₄*ρ⁴)")
    print("where ρ = ζR (ζ is the effective nuclear charge, and R is the internuclear distance).\n")
    
    # Print the values of each coefficient
    print("The coefficients in the equation are:")
    print(f"c₀ = {c0}")
    print(f"c₁ = 1/2 = {c1}")
    print(f"c₂ = 1/12 ≈ {c2:.5f}")
    print(f"c₃ = {c3}")
    print(f"c₄ = 1/240 ≈ {c4:.5f}\n")

    # Print the final, complete formula
    print("The final expression for S in terms of R and ζ is:")
    print(f"S(R, ζ) = exp(-ζ*R/2) * ({c0} + ({c1})*(ζ*R) + {c2:.5f}*(ζ*R)² + {c4:.5f}*(ζ*R)⁴)")

solve_overlap_integral()
