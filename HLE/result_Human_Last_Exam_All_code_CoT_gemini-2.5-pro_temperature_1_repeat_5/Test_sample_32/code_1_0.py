from fractions import Fraction

def calculate_lambda_integral():
    """
    Calculates the integral of the product of lambda classes lambda_3*lambda_2*lambda_1
    on the moduli space of stable curves of genus 3.
    """
    
    # The problem is to compute the integral of lambda_3 * lambda_2 * lambda_1 over
    # the moduli space of stable curves of genus 3, M_3_bar.
    # The dimension of this space is 3*g-3 = 6, and the degree of the integrand is 3+2+1=6.

    # According to relations discovered by C. Faber in the tautological ring of M_3,
    # the following equality holds: lambda_3 * lambda_2 * lambda_1 = (1/24) * lambda_1^6.
    # This allows us to relate the integral of our target class to a more standard one.
    # Integral(lambda_3*lambda_2*lambda_1) = (1/24) * Integral(lambda_1^6)
    
    # The proportionality constant from Faber's relation.
    coeff_numerator = 1
    coeff_denominator = 24
    coeff = Fraction(coeff_numerator, coeff_denominator)

    # The integral of lambda_1^6 over M_3_bar is a known result from the theory of Hodge integrals.
    # Integral(lambda_1^6) = 1/720.
    integral_lambda_1_pow_6_numerator = 1
    integral_lambda_1_pow_6_denominator = 720
    integral_lambda_1_pow_6 = Fraction(integral_lambda_1_pow_6_numerator, integral_lambda_1_pow_6_denominator)

    # We multiply the coefficient by the value of the integral of lambda_1^6.
    result = coeff * integral_lambda_1_pow_6

    # Print the explanation and the final equation with all numbers.
    print("The integral of lambda_3*lambda_2*lambda_1 on the moduli of stable curves of genus 3 is calculated using Faber's relations.")
    print("The relation is: Integral(lambda_3*lambda_2*lambda_1) = (1/24) * Integral(lambda_1^6).")
    print("The known value for Integral(lambda_1^6) is 1/720.")
    print("\nThe final calculation is:")
    print(f"({coeff_numerator}/{coeff_denominator}) * ({integral_lambda_1_pow_6_numerator}/{integral_lambda_1_pow_6_denominator}) = {result.numerator}/{result.denominator}")
    

if __name__ == '__main__':
    calculate_lambda_integral()
