# Define the operation
def operation(a, b):
    return (a + b) / b

# Calculate 5@2
result_5_2 = operation(5, 2)

# Calculate (5@2)@6
final_result = operation(result_5_2, 6)

# Print the final result
print(final_result)