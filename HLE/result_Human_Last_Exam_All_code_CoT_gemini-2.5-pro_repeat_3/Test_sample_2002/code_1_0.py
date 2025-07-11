def find_winning_t_for_p2(limit):
    """
    Finds the initial number of tokens T for which the second player has a winning strategy.
    These are the Fibonacci numbers.
    """
    winning_t_values = set()
    a, b = 0, 1
    while b <= limit:
        # The positive Fibonacci numbers are the losing positions for the first player.
        # Standard sequence: 0, 1, 1, 2, 3, 5, ...
        # The set of pile sizes T should not have duplicates.
        if b > 0:
            winning_t_values.add(b)
        a, b = b, a + b
    
    # Return a sorted list of the unique Fibonacci numbers.
    return sorted(list(winning_t_values))

# Set a limit for T to demonstrate the solution
T_LIMIT = 200

# Get the list of T values
p2_winning_t = find_winning_t_for_p2(T_LIMIT)

print(f"The values of T up to {T_LIMIT} for which the second player has a winning strategy are:")
# The problem asks to output each number in the final equation.
# We interpret this as printing the found numbers.
for t in p2_winning_t:
    print(t)