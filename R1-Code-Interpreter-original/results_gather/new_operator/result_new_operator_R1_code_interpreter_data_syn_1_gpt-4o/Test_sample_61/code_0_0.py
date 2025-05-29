# Define the operation a ~ b
def operation(a, b):
    return 3 * a + 3 * b

# Calculate 9 ~ 4
result_9_4 = operation(9, 4)

# Calculate (9 ~ 4) ~ 9
final_result = operation(result_9_4, 9)

# Print the final result
print(final_result)