def calculate_T(n):
    """
    Calculates the minimum number of trials T(n).
    The formula is T(n) = 2n - s2(n), where s2(n) is the number of 1s 
    (sum of bits) in the binary representation of n.
    """
    # Count the number of '1's in the binary representation of n
    num_set_bits = bin(n).count('1')
    
    # Calculate T(n) using the derived formula
    result = 2 * n - num_set_bits
    
    # Print the equation and the result as requested, showing each number
    print(f"T({n}) = 2 * {n} - {num_set_bits} = {result}")

# The values of n for which to calculate T(n)
n_values = [2, 3, 1234, 6712]

# Calculate and print T(n) for each value in the list
for n in n_values:
    calculate_T(n)