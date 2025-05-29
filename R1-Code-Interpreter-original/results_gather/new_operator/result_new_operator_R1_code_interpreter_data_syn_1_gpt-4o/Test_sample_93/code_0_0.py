# Define the operation a@b
def operation(a, b):
    return (a + b) / b

# Calculate 9@4
result_9_4 = operation(9, 4)

# Calculate (9@4)@3
final_result = operation(result_9_4, 3)

# Print the final result
print(final_result)