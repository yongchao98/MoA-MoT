import sympy
import math

def calculate_lambda_integral(g):
    """
    Calculates the integral of the product of lambda classes lambda_g * ... * lambda_1
    on the moduli space of stable curves of genus g.

    The formula is: |B_{2g}| / (2^(g-1) * (2g)!)
    """

    # The index for the Bernoulli number and factorial
    n = 2 * g

    # Calculate the components of the formula
    try:
        B_n = sympy.bernoulli(n)
        abs_B_n = abs(B_n)
        B_n_num, B_n_den = sympy.fraction(abs_B_n)
    except Exception as e:
        print(f"Error calculating Bernoulli number: {e}")
        return

    power_of_2 = 2**(g-1)
    
    try:
        factorial_term = math.factorial(n)
    except ValueError:
        print(f"Error: Input to factorial must be a non-negative integer.")
        return

    # Assemble the final result
    final_numerator = B_n_num
    final_denominator = B_n_den * power_of_2 * factorial_term

    # Print the breakdown of the calculation as requested
    print(f"The integral is calculated using the formula: |B_{2g}| / (2^(g-1) * (2g)!)")
    print(f"For genus g = {g}:")
    print(f"The Bernoulli number B_{n} is {B_n}.")
    print(f"The term 2^(g-1) is 2^({g-1}) = {power_of_2}.")
    print(f"The factorial (2g)! is {n}! = {factorial_term}.")
    print("\nThe final equation with the computed values is:")
    print(f"|{B_n}| / ({power_of_2} * {factorial_term}) = ({B_n_num}/{B_n_den}) / {power_of_2 * factorial_term} = {final_numerator}/{final_denominator}")

# For the specific problem of genus 3
genus = 3
calculate_lambda_integral(genus)
