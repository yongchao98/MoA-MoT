from sympy import isprime

# Possible number of true statements
possible_true_counts = [2, 3, 5, 7]  # Prime numbers

# Possible number of false statements
possible_false_counts = [4, 6]  # Composite numbers

# Total number of statements
total_statements = 7

# Check each combination of true and false counts
possible_solutions = 0
for true_count in possible_true_counts:
    false_count = total_statements - true_count
    if false_count in possible_false_counts:
        # Check the conditions
        if (true_count == 7 and false_count == 0) or (true_count == 2 and false_count == 5):
            # Statement 5 condition
            if (true_count == 7 and false_count == 0) or (true_count == 2 and false_count == 5):
                possible_solutions += 1

print(possible_solutions)