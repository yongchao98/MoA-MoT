import math

def calculate_cn_prefactor(n):
    """
    Calculates and prints the prefactor c_n for a given integer n.
    The formula for the prefactor is c_n = -(n-1) / n!

    This prefactor corresponds to the contribution of the fully f-connected
    graph (complete graph K_n) to the n-th virial coefficient B_n.

    Args:
        n (int): The number of particles/nodes in the diagram. Must be >= 2.
    """
    if not isinstance(n, int) or n < 2:
        print(f"Error: The number of nodes 'n' must be an integer greater than or equal to 2.")
        return

    # The final equation for the prefactor c_n
    final_equation_str = f"c_{n} = -({n}-1) / {n}!"
    print(f"For n = {n}:")
    print(f"The equation for the prefactor is: {final_equation_str}")

    # Calculate the numbers in the final equation
    # Numerator is -(n-1)
    numerator = -(n - 1)
    # Denominator is n!
    denominator = math.factorial(n)

    print(f"Numerator = -({n} - 1) = {numerator}")
    print(f"Denominator = {n}! = {denominator}")
    
    # Calculate the final value and print it as a fraction
    print(f"c_{n} = {numerator} / {denominator}")

    # To provide a cleaner look, we can simplify the fraction
    if numerator != 0:
        common_divisor = math.gcd(abs(numerator), denominator)
        reduced_num = numerator // common_divisor
        reduced_den = denominator // common_divisor
        if reduced_den == 1:
            print(f"As a simplified value, c_{n} = {reduced_num}")
        else:
            print(f"As a simplified fraction, c_{n} = {reduced_num}/{reduced_den}")
    else: # This case happens for n=1, but we restricted n>=2
        print(f"As a simplified value, c_{n} = 0")
        
    print("-" * 30)

if __name__ == '__main__':
    # Demonstrate the calculation for a few values of n
    calculate_cn_prefactor(2)
    calculate_cn_prefactor(3)
    calculate_cn_prefactor(4)
    calculate_cn_prefactor(5)