# Define the operation a@b
def operation(a, b):
    return (a + b) / b

# Calculate 4@9
result_4_9 = operation(4, 9)

# Calculate (4@9)@3
final_result = operation(result_4_9, 3)

# Print the final result
print(final_result)