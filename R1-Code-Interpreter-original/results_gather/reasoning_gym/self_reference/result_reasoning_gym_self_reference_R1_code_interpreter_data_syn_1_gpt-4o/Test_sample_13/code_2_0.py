# Calculate the number of possible solutions
# True statements: 4, 5, 6
# False statements: 1, 2, 3, 7

true_statements = [4, 5, 6]
false_statements = [1, 2, 3, 7]

# Check if the number of true statements is 3 and false statements is 4
num_true = len(true_statements)
num_false = len(false_statements)

# Verify the conditions
is_valid_solution = (num_true == 3) and (num_false == 4)

# Output the number of valid solutions
print(1 if is_valid_solution else 0)