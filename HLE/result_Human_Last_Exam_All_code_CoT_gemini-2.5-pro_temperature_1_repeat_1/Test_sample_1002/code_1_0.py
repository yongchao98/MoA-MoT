import sys

def calculate_limit_exponent():
    """
    This function calculates the specified limit based on the integer k.
    The problem asks for the value of lim_{m -> infinity} ln(f(m))/ln(m).
    This limit is the exponent p in the asymptotic relation f(m) = Theta(m^p).

    Based on the analysis of the extremal case (an all-ones matrix), the exponent is
    p = 1 - 1/(2k).
    """
    try:
        # Prompt the user to enter the value of k.
        k_str = input("Enter the integer k (where k >= 2): ")
        k = int(k_str)

        # Validate that k is at least 2.
        if k < 2:
            print(f"Error: The integer k must be greater than or equal to 2. You entered {k}.", file=sys.stderr)
            return

    except ValueError:
        print(f"Error: Invalid input. Please enter an integer.", file=sys.stderr)
        return
    except EOFError:
        print("Error: No input received.", file=sys.stderr)
        return

    # Calculate the numerator and denominator of the resulting fraction.
    # The formula is 1 - 1/(2k), which simplifies to (2k-1)/(2k).
    numerator = 2 * k - 1
    denominator = 2 * k

    # Print the result, showing the numbers plugged into the formula as requested.
    print(f"For k = {k}, the value of the limit is calculated as follows:")
    print(f"1 - 1 / (2 * {k}) = {numerator}/{denominator}")

if __name__ == '__main__':
    calculate_limit_exponent()
