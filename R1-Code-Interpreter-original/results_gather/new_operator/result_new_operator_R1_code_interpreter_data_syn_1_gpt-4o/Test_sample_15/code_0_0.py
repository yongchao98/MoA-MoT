# Define the operation a ~ b
def operation(a, b):
    return 3 * a + 2 * b

# Calculate 10 ~ 2
first_result = operation(10, 2)

# Calculate (10 ~ 2) ~ 2
final_result = operation(first_result, 2)

# Print the final result
print(final_result)