# Define the operation a $ b
def dollar_operation(a, b):
    return 3 * a + 3 * b

# Calculate 5 $ 9
result_5_9 = dollar_operation(5, 9)

# Calculate (5 $ 9) $ 4
final_result = dollar_operation(result_5_9, 4)

# Print the final result
print(final_result)