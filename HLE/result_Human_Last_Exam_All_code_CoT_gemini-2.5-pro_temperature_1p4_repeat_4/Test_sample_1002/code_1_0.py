import sys

def solve_limit():
    """
    This function calculates the limit lim_{m->inf} ln(f(m))/ln(m) based on k.
    """
    try:
        # Prompt the user to enter the value of k.
        k_str = input("Enter the integer k (where k >= 2): ")
        k = int(k_str)

        if k < 2:
            print("Error: k must be an integer greater than or equal to 2.", file=sys.stderr)
            return

        # The derivation shows that the limit is 1 - 1/(2k).
        # We will now print the formula and the result.
        
        # Define the components of the equation
        one_val = 1
        two_val = 2
        
        # Calculate the result
        result = 1 - 1 / (2 * k)

        # Output the formula with each number explicitly stated as requested.
        print(f"The formula for the limit is: {one_val} - {one_val} / ({two_val} * k)")
        
        # Output the final computed value
        print(f"For k = {k}, the value of the limit is: {result}")

    except ValueError:
        print("Error: Invalid input. Please enter an integer.", file=sys.stderr)
    except Exception as e:
        print(f"An unexpected error occurred: {e}", file=sys.stderr)

if __name__ == '__main__':
    solve_limit()
