# Define the operation a ~ b
def operation(a, b):
    return 2 * a + 4 * b

# Calculate 2 ~ 9
first_result = operation(2, 9)

# Calculate (2 ~ 9) ~ 9
final_result = operation(first_result, 9)

# Print the final result
print(final_result)