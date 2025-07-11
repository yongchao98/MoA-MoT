def solve_complex_integral():
    """
    Calculates the integral over contour gamma based on the given integrals
    over gamma_1 and gamma_2 using the deformation principle.
    """
    # The given integral values. In Python, the imaginary unit is represented by 'j'.
    integral_gamma1 = 3 + 4j
    integral_gamma2 = 5 + 6j

    # According to the deformation principle, the integral over gamma is the sum
    # of the integrals over gamma_1 and gamma_2 because gamma encloses both
    # singularities with the same orientation.
    result_integral = integral_gamma1 + integral_gamma2

    # For display purposes, we use 'i' as in the problem statement.
    val1_str = f"{int(integral_gamma1.real)} + {int(integral_gamma1.imag)}i"
    val2_str = f"{int(integral_gamma2.real)} + {int(integral_gamma2.imag)}i"
    
    print("The governing equation from the deformation principle is:")
    print("∫γ f dz = ∫γ₁ f dz + ∫γ₂ f dz")
    
    print("\nSubstituting the given values into the equation:")
    print(f"∫γ f dz = ({val1_str}) + ({val2_str})")
    
    # Show the addition of real and imaginary parts separately.
    real_part_1 = int(integral_gamma1.real)
    imag_part_1 = int(integral_gamma1.imag)
    real_part_2 = int(integral_gamma2.real)
    imag_part_2 = int(integral_gamma2.imag)

    print("\nAdding the real and imaginary parts:")
    print(f"∫γ f dz = ({real_part_1} + {real_part_2}) + ({imag_part_1} + {imag_part_2})i")
    
    # Display the final calculated result.
    final_real = int(result_integral.real)
    final_imag = int(result_integral.imag)
    print("\nThe final result is:")
    print(f"∫γ f dz = {final_real} + {final_imag}i")

solve_complex_integral()