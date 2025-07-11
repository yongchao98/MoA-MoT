def solve_T(n):
    """
    Calculates T(n), the minimum number of trials needed to decide if we have
    an equal number of real and fake golden bars among 2n bars.
    """
    if not isinstance(n, int) or n <= 0:
        raise ValueError("n must be a positive integer.")
    
    # Handle the special case for n=1
    if n == 1:
        return 1
    # If n is even (and n>=2)
    elif n % 2 == 0:
        return 2 * n - 1
    # If n is odd (and n>=3)
    else:
        return 2 * n - 2

# Values of n to be calculated
n_values = [2, 3, 1234, 6712]

# Calculate the results for each n
results = [solve_T(n) for n in n_values]

# Print the final answer in the required format
# "Remember in the final code you still need to output each number in the final equation!"
# The problem asks for the values separated by comma.
print(f"{results[0]},{results[1]},{results[2]},{results[3]}")