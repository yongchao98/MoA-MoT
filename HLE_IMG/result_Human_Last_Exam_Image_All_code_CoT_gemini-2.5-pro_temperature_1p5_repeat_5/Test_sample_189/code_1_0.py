def solve_complex_integral():
    """
    Calculates the integral of f over the contour gamma by summing the
    integrals over gamma_1 and gamma_2, based on the deformation of
    contours principle in complex analysis.
    """

    # The integral over gamma_1 is given as 3 + 4i.
    # In Python, the imaginary unit is represented by 'j'.
    integral_gamma1 = 3 + 4j

    # The integral over gamma_2 is given as 5 + 6i.
    integral_gamma2 = 5 + 6j

    # According to the deformation of contours principle, the integral over gamma
    # is the sum of the integrals over the contours it encloses.
    integral_gamma = integral_gamma1 + integral_gamma2

    # Print the explanation and the final equation.
    print("Based on the principle of deformation of contours, the integral over γ is the sum of the integrals over γ₁ and γ₂.")
    print("∫γ f(z)dz = ∫γ₁ f(z)dz + ∫γ₂ f(z)dz")
    # Using python's format to display the complex numbers correctly
    # and ensuring that we show the components of the sum.
    print(f"∫γ f(z)dz = ({integral_gamma1}) + ({integral_gamma2}) = {integral_gamma}")

    # To show the real and imaginary parts separately
    real_part = integral_gamma.real
    imag_part = integral_gamma.imag
    print(f"The result is {int(real_part)} + {int(imag_part)}i.")

# Execute the function
solve_complex_integral()