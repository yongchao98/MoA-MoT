import numpy as np

def solve_integral_approximation():
    """
    Calculates the coefficients for the asymptotic expansion of the integral I(epsilon)
    and prints the analytical formula.
    """

    # The denominator polynomial is f(x) = 9*x**5 + 5*x**6 + 9*x**8
    # For small x, the dominant term is c1*x**n1 where c1=9.0, n1=5.0
    # The next dominant term is c2*x**n2 where c2=5.0, n2=6.0
    n1 = 5.0
    c1_coeff = 9.0
    
    # The power of epsilon in the leading term is p1 = (n1-1)/n1
    p1 = (n1 - 1) / n1

    # The coefficient of the leading term is C1 = (pi/n1) / (c1**(1/n1) * sin(pi/n1))
    C1 = (np.pi / n1) / (c1_coeff**(1/n1) * np.sin(np.pi / n1))

    # The power of epsilon in the first correction term is p2 = (n1-2)/n1
    p2 = (n1 - 2) / n1
    
    # The coefficient of the first correction term can be derived as
    # C2 = -c2 * c1**(-(n1+1)/n1) * integral(u**(n1) / (1+u**n1)**2 du)
    # The integral evaluates to (1/n1) * Gamma((n1+1)/n1)*Gamma(1-(1/n1))
    # which simplifies using Gamma properties. A simpler form is:
    # C2 = - (c2 / c1_coeff) * ( (2*pi)/(n1*n1*sin(2*pi/n1)) ) / (c1_coeff**(1/n1))
    # where c2=5.0. An even more direct derivation yields:
    C2 = - (2 * np.pi) / (n1 * c1_coeff**((n1+2)/n1) * np.sin(2*np.pi/n1))
    
    # In our case c2=5.0 appears in the numerator. Let's rederive it to be sure.
    # C2_term_coefficient = -c2_coeff * c1_coeff**(-(n1+2)/n1)
    # Integral part = integral u**(n1+1) / (1+u**n1)**2 du
    # let t = u**n1 -> integral = (1/n1)*integral t**(2/n1)/(1+t)**2 dt = (1/n1)B(2/n1+1, 1-2/n1)
    # The integral term is (1/5) * Gamma(7/5)Gamma(3/5) = (1/5) * (2/5)Gamma(2/5)Gamma(3/5) = (2/25)*pi/sin(2pi/5)
    # So C2 = -5 * 9**(-7/5) * (2/25)*pi/sin(2pi/5) = -(2*pi)/(5*9**(7/5)*sin(2*pi/5))
    c2_coeff_from_problem = 5.0
    C2_final = - (c2_coeff_from_problem/c1_coeff) * (2 * np.pi) / (n1 * c1_coeff**(2./n1) * np.sin(2 * np.pi/n1))

    print("The analytical approximation for the integral I(epsilon) for small epsilon is:")
    print(f"I(eps) ≈ C1 * eps**(-p1) + C2 * eps**(-p2)")
    print("where:")
    print(f"p1 = {p1:.1f}")
    print(f"C1 = {C1:.4f}")
    print(f"p2 = {p2:.1f}")
    print(f"C2 = {C2_final:.4f}")
    print("\nFinal formula:")
    print(f"I(eps) ≈ {C1:.4f} * eps**(-{p1:.1f}) {C2_final:+.4f} * eps**(-{p2:.1f})")


solve_integral_approximation()