import math

def calculate_function_field_mertens_limit(q):
    """
    Calculates the value of the limit for a given characteristic q.
    
    The formula is (e^-gamma * log(q)) / (q - 1).
    """
    
    # The Euler-Mascheroni constant
    euler_gamma = 0.57721566490153286060651209008240243104215933593992

    # Calculate the components of the formula
    e_to_minus_gamma = math.exp(-euler_gamma)
    log_q = math.log(q)
    q_minus_1 = q - 1

    # Calculate the final value
    result = (e_to_minus_gamma * log_q) / q_minus_1

    # Print the components and the final result
    print(f"For q = {q}:")
    print(f"  Euler-Mascheroni constant (gamma): {euler_gamma}")
    print(f"  The term e^(-gamma): {e_to_minus_gamma}")
    print(f"  The term log(q): {log_q}")
    print(f"  The term q - 1: {q_minus_1}")
    print("\nFinal equation: (e^(-gamma) * log(q)) / (q - 1)")
    print(f"Result: ({e_to_minus_gamma} * {log_q}) / {q_minus_1}")
    print(f"\nThe value of the limit is: {result}")
    
    return result

if __name__ == "__main__":
    # The characteristic q must be a prime power, q > 1.
    # As an example, we use q=5.
    q_example = 5
    final_value = calculate_function_field_mertens_limit(q_example)
    # The final answer is the numerical value printed above.
    # To conform to the requested format, we output the value in the specified tag.
    # The output format requires a single value.
    # print(f"<<<{final_value}>>>") # This would be the final line if required.
