import sys
import math

def calculate_prefactor_cn(n_str):
    """
    This script calculates the system-independent prefactor c_n for the
    fully f-connected Ree-Hoover diagram contribution to the n-th virial coefficient.

    The prefactor is given by the formula: c_n = -(n-1) / n!
    """
    try:
        n = int(n_str)
        if n < 2:
            print("Error: Virial coefficients are defined for n >= 2.", file=sys.stderr)
            sys.exit(1)
    except ValueError:
        print(f"Error: Invalid input '{n_str}'. Please provide an integer.", file=sys.stderr)
        sys.exit(1)

    # Calculate the components of the formula
    numerator = n - 1
    denominator = math.factorial(n)
    
    # Calculate the final result
    result = -numerator / denominator

    # Print the explanation and the step-by-step calculation as requested
    print(f"To find the prefactor c_n for n = {n}, we use the formula:")
    print(f"c_n = -(n - 1) / n!")
    print("\nSubstituting the value of n:")
    print(f"c_{n} = -({n} - 1) / {n}!")
    print(f"c_{n} = -{numerator} / {denominator}")
    print(f"\nThe calculated value is:")
    print(f"c_{n} = {result}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(f"Usage: python {sys.argv[0]} <n>", file=sys.stderr)
        print("Please provide an integer n (>= 2) as a command-line argument.", file=sys.stderr)
        sys.exit(1)
    
    calculate_prefactor_cn(sys.argv[1])