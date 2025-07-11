import math

def calculate_u_r(n):
    """
    Calculates the minimal order of the Picard-Fuchs differential equation u_r(n).
    The formula is u_r(n) = 2 * floor((n-1)/2).
    """
    # In Python, integer division // is equivalent to math.floor() for positive results.
    return 2 * ((n - 1) // 2)

def solve_task():
    """
    Solves the main task by calculating u_r(n) for n from 3 to 12
    and printing the results.
    """
    # The range of n for which to calculate u_r(n)
    n_values = range(3, 13)

    # Calculate the values
    u_r_values = [calculate_u_r(n) for n in n_values]

    # Create the output string representation of the set of equations
    # e.g., "{u_r(3) = 2, u_r(4) = 2, ...}"
    equations = [f"u_r({n}) = {val}" for n, val in zip(n_values, u_r_values)]
    output_string = "{" + ", ".join(equations) + "}"

    # Print the final result, showing each number in the final equation set
    print(f"The set of values for {{u_r(3), u_r(4), ..., u_r(12)}} is:")
    print(output_string)

solve_task()