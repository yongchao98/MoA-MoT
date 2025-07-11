import math
import sys

def calculate_and_print_cn(n):
    """
    Calculates and prints the prefactor c_n for a given integer n >= 2.
    
    The prefactor c_n appears in the virial expansion and is associated
    with the fully f-connected Ree-Hoover diagram. The derived formula is:
    c_n = - (n - 1) / n!
    
    Args:
        n (int): The number of particles/nodes, must be >= 2.
    """
    # --- Input validation ---
    if not isinstance(n, int) or n < 2:
        print(f"Error: n must be an integer greater than or equal to 2, but got {n}.", file=sys.stderr)
        return

    # --- Calculation Steps ---
    # 1. Calculate the numerator: -(n - 1)
    numerator = -(n - 1)
    
    # 2. Calculate the denominator: n!
    denominator = math.factorial(n)

    # 3. Calculate the final value of c_n
    cn_value = numerator / denominator

    # --- Output the results ---
    print(f"Calculation of the prefactor c_n for n = {n}:")
    print("-" * 40)
    
    # Final Equation is c_n = -(n - 1) / n!
    print("The derived formula for c_n is: -(n - 1) / n!")
    
    # Show the formula with the specific number for n
    print(f"Substituting n = {n}:")
    print(f"c_{n} = -({n} - 1) / {n}!")
    
    # Show the values of the numerator and denominator
    print(f"c_{n} = {numerator} / {denominator}")
    
    # Show the final numerical result
    print(f"c_{n} = {cn_value}")
    print("-" * 40)


def main():
    """
    Main function to demonstrate the calculation of c_n.
    """
    # The problem asks to determine the general prefactor c_n.
    # The code calculates and displays this value for a specific choice of n.
    # We choose n=4 as an example.
    
    # You can change this value to calculate c_n for other n.
    n_example = 4
    calculate_and_print_cn(n_example)
    
    # Another example for n=5
    # n_example = 5
    # calculate_and_print_cn(n_example)

if __name__ == "__main__":
    main()