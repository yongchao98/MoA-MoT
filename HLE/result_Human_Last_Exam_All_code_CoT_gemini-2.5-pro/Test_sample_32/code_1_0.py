from fractions import Fraction

def compute_integral():
    """
    Computes the integral of lambda_3*lambda_2*lambda_1 on M_3.
    
    The calculation is based on the following facts from the intersection theory
    on the moduli space of curves:
    1. The relation in the tautological ring R^6(M_3): 
       lambda_3 * lambda_2 * lambda_1 = (5/12) * kappa_3^2
    2. The known integral of the basis element:
       integral(kappa_3^2) over M_3 = 1/720
    """
    
    # The proportionality constant in the relation lambda_3*lambda_2*lambda_1 = c * kappa_3^2
    c = Fraction(5, 12)
    
    # The known value of the integral of kappa_3^2
    integral_kappa_3_squared = Fraction(1, 720)
    
    # The integral is the product of the constant and the known integral
    result = c * integral_kappa_3_squared
    
    # Print the equation with all the numbers
    print(f"The integral is calculated by the formula: integral(lambda_3*lambda_2*lambda_1) = c * integral(kappa_3^2)")
    print(f"c = {c.numerator}/{c.denominator}")
    print(f"integral(kappa_3^2) = {integral_kappa_3_squared.numerator}/{integral_kappa_3_squared.denominator}")
    print(f"Result: ({c.numerator}/{c.denominator}) * ({integral_kappa_3_squared.numerator}/{integral_kappa_3_squared.denominator}) = {result.numerator}/{result.denominator}")
    print("\nFinal answer:")
    print(f"{result.numerator}/{result.denominator}")

compute_integral()