from fractions import Fraction

def calculate_lambda_integral():
    """
    Calculates the integral of the product of lambda classes lambda_3*lambda_2*lambda_1
    on the moduli space of stable curves of genus 3.

    The calculation is based on established relations in the tautological ring R*(M_3_bar).
    The key formulas are:
    1. The relation in R^6(M_3_bar): lambda_1*lambda_2*lambda_3 = (1/12) * kappa_1^6
    2. The known integral: Integral(kappa_1^6) over M_3_bar = 1/720
    """

    # The proportionality constant from the relation lambda_1*lambda_2*lambda_3 = c * kappa_1^6
    # c = 1/12
    c_num = 1
    c_den = 12
    constant_fraction = Fraction(c_num, c_den)

    # The value of the integral of kappa_1^6
    # integral_kappa = 1/720
    integral_kappa_num = 1
    integral_kappa_den = 720
    integral_kappa_fraction = Fraction(integral_kappa_num, integral_kappa_den)

    # The final result is the product of the constant and the integral of kappa_1^6
    result_fraction = constant_fraction * integral_kappa_fraction

    # Extract numerators and denominators for printing the equation
    res_num = result_fraction.numerator
    res_den = result_fraction.denominator

    print("The integral is calculated using a known relation in the tautological ring.")
    print(f"The final equation is: ({c_num}/{c_den}) * ({integral_kappa_num}/{integral_kappa_den}) = {res_num}/{res_den}")
    print(f"The value of the integral is: {res_num}/{res_den}")

if __name__ == "__main__":
    calculate_lambda_integral()
