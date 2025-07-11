def solve_integral():
    """
    Calculates the integral over contour gamma based on the integrals over
    gamma_1 and gamma_2.
    """
    # Given integral values
    integral_gamma1 = complex(3, 4)  # ∫γ₁ f = 3 + 4i
    integral_gamma2 = complex(5, 6)  # ∫γ₂ f = 5 + 6i

    # The contour γ goes around z₂ in the same direction as γ₂ (counter-clockwise)
    # but around z₁ in the opposite direction of γ₁ (clockwise).
    # So, ∫γ f = ∫γ₂ f - ∫γ₁ f
    result = integral_gamma2 - integral_gamma1

    # Extract integer parts for clean printing
    g2_real = int(integral_gamma2.real)
    g2_imag = int(integral_gamma2.imag)
    g1_real = int(integral_gamma1.real)
    g1_imag = int(integral_gamma1.imag)
    res_real = int(result.real)
    res_imag = int(result.imag)
    
    # Print the full calculation as requested.
    print(f"The relationship between the integrals is: ∫γ f = ∫γ₂ f - ∫γ₁ f")
    print(f"Substituting the given values, we get:")
    # Using python's f-string formatting to display the equation clearly.
    # We display each number in the final equation.
    print(f"({g2_real} + {g2_imag}i) - ({g1_real} + {g1_imag}i) = {res_real} + {res_imag}i")
    
    # Final answer for verification, matching the last print statement.
    final_answer = f"{res_real} + {res_imag}i"
    # This part of the code won't print but holds the string for the <<<...>>> answer.

solve_integral()