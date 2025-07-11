import math

# Set the number of symbols, n. You can change this value to solve for different n.
n = 4

def telephone_number(num):
    """
    Calculates the n-th telephone number (involution number), T_n,
    using the recurrence relation T_n = T_{n-1} + (n-1)*T_{n-2}.
    """
    if not isinstance(num, int) or num < 0:
        raise ValueError("Input must be a non-negative integer.")
    
    # Base cases
    if num == 0:
        return 1
    if num == 1:
        return 1
    
    # Use a list to store values for T_0, T_1, ..., T_n for dynamic programming
    t = [0] * (num + 1)
    t[0] = 1
    t[1] = 1
    
    # Iteratively compute T_i up to T_n
    for i in range(2, num + 1):
        t[i] = t[i-1] + (i - 1) * t[i-2]
        
    return t[num]

# --- Main execution ---
if isinstance(n, int) and n >= 0:
    # Step 1: Calculate n!
    n_factorial = math.factorial(n)
    
    # Step 2: Calculate the n-th telephone number, T_n
    tn_value = telephone_number(n)
    
    # Step 3: Calculate the total number of cool strings
    total_strings = n_factorial * tn_value
    
    print(f"For n = {n} symbols, the number of cool strings of maximal length 3n is calculated.")
    print("The formula is: N(n) = n! * T_n")
    print("Here, n! is the factorial of n, and T_n is the n-th telephone number (or involution number).")
    print("-" * 20)
    print(f"For n = {n}:")
    print(f"n! = {n_factorial}")
    print(f"T_{n} = {tn_value}")
    print(f"The final calculation is: {n_factorial} * {tn_value} = {total_strings}")
else:
    print("Error: The number of symbols 'n' must be a non-negative integer.")
