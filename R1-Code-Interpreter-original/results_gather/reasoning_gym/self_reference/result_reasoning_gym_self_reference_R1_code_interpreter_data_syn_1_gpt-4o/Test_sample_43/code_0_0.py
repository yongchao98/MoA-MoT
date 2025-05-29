# Define the statements
def check_statements():
    # Possible number of true statements
    possible_true_counts = [4, 5, 6, 7]
    # Possible number of false statements
    possible_false_counts = [0, 1, 2, 3, 4, 5]
    # Prime numbers for true statements
    prime_true_counts = [2, 3, 5, 7]
    # Composite numbers for false statements
    composite_false_counts = [4, 6]

    # Check each possible number of true statements
    for T in possible_true_counts:
        F = 7 - T
        # Check if the scenario satisfies all conditions
        if (T >= 4 and F <= 5 and
            (T == 0 or F == 3) and
            (T in prime_true_counts) and
            (F in composite_false_counts)):
            return 1  # Found a valid solution

    return 0  # No valid solution found

# Execute the function and print the result
print(check_statements())