import math

def calculate_and_print_cn(n):
    """
    Calculates the prefactor c_n for a given integer n >= 2.

    The prefactor c_n relates the nth virial coefficient B_n to the integral
    Lambda_n over the fully f-connected diagram (a complete graph K_n).
    The formula is c_n = -(n-1) / n!

    Args:
        n (int): The number of particles (nodes in the diagram).
    """
    if not isinstance(n, int) or n < 2:
        print(f"Error: Input n must be an integer greater than or equal to 2.")
        return

    # Calculate the numerator and denominator of the fraction
    numerator = n - 1
    denominator = math.factorial(n)

    # Calculate the final value of c_n
    cn_value = -numerator / denominator

    # Print the results, showing each number in the final equation
    print(f"For n = {n}:")
    print(f"  The formula for the prefactor is: c_n = - (n - 1) / n!")
    print(f"  Substituting the value of n:")
    print(f"  c_{n} = - ({n} - 1) / {n}!")
    print(f"  c_{n} = - {numerator} / {denominator}")
    print(f"  The numerical value is: c_{n} = {cn_value}\n")

if __name__ == '__main__':
    print("This script determines the system-independent prefactor c_n.")
    print("The prefactor arises in the decomposition of the virial coefficient B_n,")
    print("isolating the contribution from the fully f-connected diagram.\n")

    # Demonstrate the calculation for a few values of n
    for n_value in range(2, 6):
        calculate_and_print_cn(n_value)
