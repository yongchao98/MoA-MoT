import math

def calculate_flagstones():
    """
    Calculates the number of flagstones for the Theatre Square problem
    with the inputs from Question 4 and prints the final equation.
    """
    # Define the input values from the problem
    n = 4 * 10**9
    m = 4 * 10**9
    a = 1

    # In Python, integers have arbitrary precision, so overflow is not a concern.
    # We calculate the number of stones needed along each dimension using ceiling division.
    flagstones_along_n = math.ceil(n / a)
    flagstones_along_m = math.ceil(m / a)

    # The total number of flagstones is the product of the two dimensions.
    total_flagstones = flagstones_along_n * flagstones_along_m

    # Print the final equation with each number, formatted as integers.
    print(f"{int(flagstones_along_n)} * {int(flagstones_along_m)} = {int(total_flagstones)}")

calculate_flagstones()