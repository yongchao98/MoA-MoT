# Define the operation
def operation(a, b):
    return 2 * a + 5 * b

# Calculate 2@4
first_result = operation(2, 4)

# Calculate (2@4)@4
final_result = operation(first_result, 4)

# Print the final result
print(final_result)