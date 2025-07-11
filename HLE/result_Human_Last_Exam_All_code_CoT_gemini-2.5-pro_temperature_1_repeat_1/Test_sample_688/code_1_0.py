import math

def calculate_and_print_cn(n):
    """
    Calculates the system-independent prefactor c_n for the fully
    f-connected Ree-Hoover diagram contribution to the n-th virial coefficient.
    
    The formula is c_n = -(n-1) / n!
    """
    # The number of nodes n must be 2 or greater for virial coefficients.
    if not isinstance(n, int) or n < 2:
        print(f"Invalid input n={n}. n must be an integer greater than or equal to 2.")
        return

    # Calculate the numerator and denominator of the prefactor c_n
    numerator = -(n - 1)
    denominator = math.factorial(n)

    # Simplify the fraction by finding the greatest common divisor (GCD)
    common_divisor = math.gcd(abs(numerator), denominator)
    simplified_numerator = numerator // common_divisor
    simplified_denominator = denominator // common_divisor

    # Print the results, showing the numbers in the equation
    print(f"For n = {n}:")
    print(f"c_{n} = -({n} - 1) / {n}!")
    print(f"c_{n} = {numerator} / {denominator}")
    
    # Print the simplified fraction
    if simplified_denominator == 1:
        print(f"Simplified: c_{n} = {simplified_numerator}")
    else:
        print(f"Simplified: c_{n} = {simplified_numerator}/{simplified_denominator}")
    print("-" * 25)

if __name__ == "__main__":
    print("This script calculates the prefactor c_n for the n-th virial coefficient B_n.")
    print("This prefactor corresponds to the fully f-connected diagram contribution.\n")
    
    # Calculate and display c_n for n from 2 to 6
    for n_val in range(2, 7):
        calculate_and_print_cn(n_val)