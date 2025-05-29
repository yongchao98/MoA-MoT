# Define the operation a & b
def operation(a, b):
    return 2 * a * b

# Calculate 7 & 5
result_7_5 = operation(7, 5)

# Calculate (7 & 5) & 9
final_result = operation(result_7_5, 9)

# Print the final result
print(final_result)