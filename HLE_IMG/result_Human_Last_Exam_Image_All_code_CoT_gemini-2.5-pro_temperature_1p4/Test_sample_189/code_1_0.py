import cmath

def solve_integral():
    """
    Calculates the integral over contour gamma based on the given integrals
    over gamma_1 and gamma_2.
    """
    # Define the complex numbers for the given integrals
    integral_gamma1 = complex(3, 4)
    integral_gamma2 = complex(5, 6)

    # According to the deformation principle in complex analysis, the integral
    # over the contour γ is the sum of the integrals over γ_1 and γ_2.
    # ∫_γ f = ∫_γ_1 f + ∫_γ_2 f
    integral_gamma = integral_gamma1 + integral_gamma2

    # Print the equation with the given values
    print("Based on the principle of deformation of contours:")
    print("∫_γ f = ∫_γ_1 f + ∫_γ_2 f")
    print(f"∫_γ f = ({int(integral_gamma1.real)} + {int(integral_gamma1.imag)}i) + ({int(integral_gamma2.real)} + {int(integral_gamma2.imag)}i)")
    
    # Print the final result
    print(f"∫_γ f = {int(integral_gamma.real)} + {int(integral_gamma.imag)}i")

solve_integral()