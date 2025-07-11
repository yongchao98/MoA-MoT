import math

def solve_theatre_square():
    """
    Solves the Theatre Square problem for the inputs given in question 4.
    """
    # Input values from question 4
    n = 4 * 10**9
    m = 4 * 10**9
    a = 1

    # To find the number of flagstones needed to cover a length, we need to find
    # the ceiling of the division. For positive integers, ceiling division
    # ceil(x / y) can be calculated using integer arithmetic as (x + y - 1) // y.
    
    # Calculate the number of flagstones needed for the 'n' dimension.
    num_a_n = (n + a - 1) // a
    
    # Calculate the number of flagstones needed for the 'm' dimension.
    num_a_m = (m + a - 1) // a

    # The total number of flagstones is the product of the two dimensions.
    total_flagstones = num_a_n * num_a_m

    # As requested, printing each number in the final equation.
    print(f"To cover dimension n ({n}), we need {num_a_n} flagstones.")
    print(f"To cover dimension m ({m}), we need {num_a_m} flagstones.")
    print("The final equation to calculate the total is:")
    print(f"{num_a_n} * {num_a_m} = {total_flagstones}")

solve_theatre_square()