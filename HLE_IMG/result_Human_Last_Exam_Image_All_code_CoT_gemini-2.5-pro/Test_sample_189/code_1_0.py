import cmath

def solve_integral():
    """
    Calculates the integral over contour gamma by summing the integrals
    over gamma1 and gamma2.
    """
    # Given values for the integrals over gamma_1 and gamma_2
    integral_gamma1 = complex(3, 4)
    integral_gamma2 = complex(5, 6)

    # According to the principle of additivity of contour integrals for the given
    # configuration, the integral over γ is the sum of the integrals over γ₁ and γ₂.
    # ∫_γ f(z)dz = ∫_γ₁ f(z)dz + ∫_γ₂ f(z)dz
    integral_gamma = integral_gamma1 + integral_gamma2

    # Print the equation with all the numbers.
    # We use int() to display the numbers without decimals for a cleaner look.
    print(f"The integral over γ is the sum of the integrals over γ₁ and γ₂.")
    print(f"The calculation is as follows:")
    print(f"({int(integral_gamma1.real)} + {int(integral_gamma1.imag)}i) + ({int(integral_gamma2.real)} + {int(integral_gamma2.imag)}i) = {int(integral_gamma.real)} + {int(integral_gamma.imag)}i")

solve_integral()