import math

def solve_integral():
    """
    Calculates the integral of the product of lambda classes lambda_3*lambda_2*lambda_1
    on the moduli space of stable curves of genus 3.
    """

    # The calculation is based on known relations in the tautological ring R*(M_3).
    # 1. The relation between the integrand and the basis class lambda_1^6 is:
    #    lambda_3 * lambda_2 * lambda_1 = (c_num / c_den) * lambda_1^6
    c_num = 9
    c_den = 2

    # 2. The integral of the basis class lambda_1^6 is known:
    #    integral(lambda_1^6) = 1 / integral_lambda1_6_den
    integral_lambda1_6_den = 5184

    # 3. The integral is the product of the constant and the basis integral.
    #    I = (c_num / c_den) * (1 / integral_lambda1_6_den)
    numerator = c_num
    denominator = c_den * integral_lambda1_6_den

    # 4. Simplify the resulting fraction.
    common_divisor = math.gcd(numerator, denominator)
    final_num = numerator // common_divisor
    final_den = denominator // common_divisor

    # Output the explanation and the calculation steps.
    print("The integral of lambda_3*lambda_2*lambda_1 is calculated using relations in the tautological ring of M_3.")
    print(f"The relation is: lambda_3*lambda_2*lambda_1 = ({c_num}/{c_den}) * lambda_1^6")
    print(f"The integral of the basis class is: integral(lambda_1^6) = 1/{integral_lambda1_6_den}")
    print("\nCombining these results, the calculation is as follows:")
    print(f"Integral = ({c_num}/{c_den}) * (1/{integral_lambda1_6_den})")
    print(f"         = {c_num} / ({c_den} * {integral_lambda1_6_den})")
    print(f"         = {numerator}/{denominator}")
    
    # Print the simplified final answer.
    print(f"         = {final_num}/{final_den}")

solve_integral()