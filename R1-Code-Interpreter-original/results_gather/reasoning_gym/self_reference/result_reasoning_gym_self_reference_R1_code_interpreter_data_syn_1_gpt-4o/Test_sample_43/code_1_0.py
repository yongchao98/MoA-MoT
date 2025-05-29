# Define the statements
def check_statements():
    # Prime numbers for true statements
    prime_true_counts = [2, 3, 5, 7]
    # Composite numbers for false statements
    composite_false_counts = [4, 6]

    # Check each possible number of true statements
    for T in range(8):  # T can be from 0 to 7
        F = 7 - T
        # Check if the scenario satisfies all conditions
        if (T >= 4 and F <= 5 and
            F == 3 and
            T in prime_true_counts and
            F in composite_false_counts):
            return 1  # Found a valid solution

    return 0  # No valid solution found

# Execute the function and print the result
print(check_statements())