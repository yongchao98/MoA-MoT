import math

def calculate_T(n):
    """
    Calculates T(n) based on the formula: T(n) = n + p2 - 1,
    where p2 is the largest power of 2 less than or equal to n.
    """
    if n <= 0:
        # T(0) is not part of the problem, but would be 0.
        return 0
    
    # Find the largest power of 2 less than or equal to n.
    # This can be calculated as 2 to the power of the floor of log base 2 of n.
    # An efficient integer-based method is to use bit manipulation.
    largest_power_of_2 = 1 << (n.bit_length() - 1)
    
    # Apply the derived formula
    result = n + largest_power_of_2 - 1
    
    # Output the calculation as requested
    print(f"T({n}) = {n} + {largest_power_of_2} - 1 = {result}")
    
    return result

def solve_puzzle():
    """
    Solves the puzzle for the given values of n and prints the results.
    """
    values_of_n = [2, 3, 1234, 6712]
    results = [calculate_T(n) for n in values_of_n]
    
    # The final answer is the sequence of results, separated by a comma.
    final_answer_string = ",".join(map(str, results))
    print(f"\nThe comma-separated values for T(2), T(3), T(1234), T(6712) are: {final_answer_string}")
    # Return the final string for the required output format.
    return final_answer_string

# Execute the solver
final_answer = solve_puzzle()

# The final answer in the requested format
# print(f"<<<{final_answer}>>>")