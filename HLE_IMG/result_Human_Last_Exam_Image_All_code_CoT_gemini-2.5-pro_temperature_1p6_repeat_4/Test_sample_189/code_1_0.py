def solve_complex_integral():
    """
    Calculates the integral over contour gamma based on the integrals over
    gamma_1 and gamma_2 using the principle of contour deformation.
    """
    # Given integrals as complex numbers
    integral_gamma1 = complex(3, 4)
    integral_gamma2 = complex(5, 6)

    # The integral over gamma is the sum of the integrals over gamma1 and gamma2
    integral_gamma = integral_gamma1 + integral_gamma2

    # Prepare numbers for the final equation output
    real1 = integral_gamma1.real
    imag1 = integral_gamma1.imag
    real2 = integral_gamma2.real
    imag2 = integral_gamma2.imag
    final_real = integral_gamma.real
    final_imag = integral_gamma.imag

    # Print the equation with all the numbers
    print(f"Based on the deformation principle, the calculation is as follows:")
    print(f"({int(real1)} + {int(imag1)}i) + ({int(real2)} + {int(imag2)}i) = {int(final_real)} + {int(final_imag)}i")

solve_complex_integral()