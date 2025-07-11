def solve_complex_integral():
    """
    Calculates the integral over a composite contour by summing the integrals
    over its component contours based on the principle of deformation.
    """
    # Define the complex numbers representing the given integrals
    integral_gamma1 = 3 + 4j
    integral_gamma2 = 5 + 6j

    # According to the principle of deformation of contours, the integral over γ
    # is the sum of the integrals over γ1 and γ2.
    integral_gamma = integral_gamma1 + integral_gamma2

    # Extract real and imaginary parts for clear printing
    real1 = int(integral_gamma1.real)
    imag1 = int(integral_gamma1.imag)
    real2 = int(integral_gamma2.real)
    imag2 = int(integral_gamma2.imag)
    real_result = int(integral_gamma.real)
    imag_result = int(integral_gamma.imag)

    # Print the equation and the step-by-step calculation
    print("The integral over the contour γ is the sum of the integrals over γ₁ and γ₂.")
    print("This is based on the deformation principle in complex analysis.")
    print(f"∫γ f = ∫γ₁ f + ∫γ₂ f")
    print(f"∫γ f = ({real1} + {imag1}i) + ({real2} + {imag2}i)")
    print(f"∫γ f = ({real1} + {real2}) + ({imag1} + {imag2})i")
    print(f"∫γ f = {real_result} + {imag_result}i")

solve_complex_integral()