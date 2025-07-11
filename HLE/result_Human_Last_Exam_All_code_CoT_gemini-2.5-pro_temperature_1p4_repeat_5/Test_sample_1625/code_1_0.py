import math

def demonstrate_fourier_series():
    """
    This function demonstrates the structure of a Fourier series, which is used
    for poloidal dependence in toroidal systems.

    We will represent a hypothetical physical quantity F(theta) as a sum of
    a few Fourier components (poloidal harmonics).
    F(theta) = 1.2 + 0.5 * cos(theta) - 0.3 * sin(theta) + 0.1 * cos(2*theta)
    """

    # Coefficients for the Fourier series
    # A0/2 is the constant term
    a0_div_2 = 1.2
    # An are the cosine coefficients for n > 0
    a_coeffs = {1: 0.5, 2: 0.1}
    # Bn are the sine coefficients for n > 0
    b_coeffs = {1: -0.3}

    print("In toroidal systems, quantities are periodic in the poloidal angle (theta).")
    print("A Fourier series is used to represent this periodic dependence.")
    print("The general form is: F(theta) = A0/2 + Sum[ An*cos(n*theta) + Bn*sin(n*theta) ]")
    print("\nAn example equation is constructed below:")

    equation_parts = []
    
    # The constant term (n=0)
    # This can also be written as a0_div_2 * cos(0*theta) since cos(0)=1
    equation_parts.append(f"{a0_div_2:.2f} * cos(0*theta)")

    # The n=1 terms
    if 1 in a_coeffs:
        sign = "+" if a_coeffs[1] >= 0 else "-"
        equation_parts.append(f"{sign} {abs(a_coeffs[1]):.2f} * cos(1*theta)")
    if 1 in b_coeffs:
        sign = "+" if b_coeffs[1] >= 0 else "-"
        equation_parts.append(f"{sign} {abs(b_coeffs[1]):.2f} * sin(1*theta)")
        
    # The n=2 terms
    if 2 in a_coeffs:
        sign = "+" if a_coeffs[2] >= 0 else "-"
        equation_parts.append(f"{sign} {abs(a_coeffs[2]):.2f} * cos(2*theta)")
    if 2 in b_coeffs:
        sign = "+" if b_coeffs[2] >= 0 else "-"
        equation_parts.append(f"{sign} {abs(b_coeffs[2]):.2f} * sin(2*theta)")

    final_equation = "F(theta) = " + " ".join(equation_parts)
    print(final_equation)

if __name__ == '__main__':
    demonstrate_fourier_series()