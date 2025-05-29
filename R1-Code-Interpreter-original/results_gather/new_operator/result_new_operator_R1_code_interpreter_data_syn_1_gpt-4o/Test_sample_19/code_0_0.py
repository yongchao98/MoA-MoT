# Define the operation
def operation(a, b):
    return 2 * a + 2 * b

# Calculate 9@4
result_9_4 = operation(9, 4)

# Calculate (9@4)@4
final_result = operation(result_9_4, 4)

# Print the final result
print(final_result)