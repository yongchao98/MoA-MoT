import sys

def calculate_sum_approximation(n_str):
    """
    Calculates the approximation for the sum based on the derived formula.
    The formula is S_n approx n^2/2 + 1/120 + 1/(252*n).
    """
    try:
        n = int(n_str)
        if n <= 0:
            print("Error: Please provide a positive integer for n.")
            return
    except ValueError:
        print("Error: Invalid input. Please provide a positive integer for n.")
        return

    # Calculate the three terms of the approximation formula
    term1 = n**2 / 2.0
    term2 = 1.0 / 120.0
    term3 = 1.0 / (252.0 * n)
    
    # Calculate the final result
    result = term1 + term2 + term3

    # Output the result, showing each term explicitly as requested
    print(f"For n = {n}, the approximation is given by the formula:")
    print(f"Approximation = n^2/2 + 1/120 + 1/(252*n)")
    print(f"\nBreaking down the calculation:")
    # Using format specifiers to control the number of decimal places for clarity
    print(f"Term 1 (n^2/2):       {term1:.10f}")
    print(f"Term 2 (1/120):       {term2:.10f}")
    print(f"Term 3 (1/(252*n)):   {term3:.10f}")
    print("---------------------------------------")
    print(f"Final approximated sum: {result:.10f}")


def main():
    """
    Main function to handle command-line arguments.
    """
    # Default value for demonstration if no argument is provided
    default_n = "100"
    
    if len(sys.argv) > 1:
        n_value = sys.argv[1]
    else:
        print(f"No value for n provided. Using default n = {default_n}.")
        n_value = default_n

    calculate_sum_approximation(n_value)

if __name__ == "__main__":
    main()
