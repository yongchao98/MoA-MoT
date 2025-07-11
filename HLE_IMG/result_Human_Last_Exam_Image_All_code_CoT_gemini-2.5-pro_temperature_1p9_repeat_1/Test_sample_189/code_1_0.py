def solve_complex_integral():
    """
    Calculates the complex integral over contour gamma by summing the
    integrals over gamma_1 and gamma_2, based on the principle of
    deformation of contours.
    """
    # Define the given values of the integrals as complex numbers.
    # In Python, the imaginary unit is represented by 'j'.
    integral_gamma_1 = 3 + 4j
    integral_gamma_2 = 5 + 6j

    # The integral over γ is the sum of the integrals over γ₁ and γ₂.
    integral_gamma = integral_gamma_1 + integral_gamma_2

    # Extract the components of each complex number for clear printing.
    val1_real = int(integral_gamma_1.real)
    val1_imag = int(integral_gamma_1.imag)

    val2_real = int(integral_gamma_2.real)
    val2_imag = int(integral_gamma_2.imag)

    result_real = int(integral_gamma.real)
    result_imag = int(integral_gamma.imag)

    # Print the equation, showing each number involved in the calculation.
    # We use 'i' for the imaginary unit in the output for standard mathematical notation.
    print(f"The integral over γ is the sum of the integrals over γ₁ and γ₂:")
    print(f"∫γ f = ({val1_real} + {val1_imag}i) + ({val2_real} + {val2_imag}i) = {result_real} + {result_imag}i")

solve_complex_integral()