import sys

def calculate_limit_value(k):
    """
    Calculates and prints the derived limit value for a given k.

    The problem asks for the value of lim_{m->inf} [ln(f(m)) / ln(m)].
    Our derivation shows this limit is equal to k / (k + 1).
    This function computes and displays this result.
    """
    if not isinstance(k, int) or k < 2:
        print("Error: Input k must be an integer greater than or equal to 2.", file=sys.stderr)
        return

    numerator = k
    denominator = k + 1
    result = numerator / denominator

    print(f"For k = {k}, the limit is given by the expression: k / (k + 1)")
    print("Substituting the value of k into the equation:")
    # As requested, printing each number in the final equation
    print(f"Numerator = k = {numerator}")
    print(f"Denominator = k + 1 = {k} + 1 = {denominator}")
    print(f"Final expression: {numerator} / {denominator}")
    print(f"The numerical value of the limit is: {result}")

def main():
    """
    Main function to handle script execution.
    It parses the command-line argument for k and calls the calculation function.
    """
    if len(sys.argv) > 1:
        try:
            k_value = int(sys.argv[1])
            calculate_limit_value(k_value)
        except ValueError:
            print(f"Error: Invalid input '{sys.argv[1]}'. Please provide an integer k >= 2.", file=sys.stderr)
    else:
        print("This script computes the value of k / (k + 1) for a given integer k >= 2.")
        print("Please provide k as a command-line argument.")
        print("For example: python your_script_name.py 3")
        print("\n--- Running an example with k = 2 ---")
        calculate_limit_value(2)

if __name__ == "__main__":
    main()
