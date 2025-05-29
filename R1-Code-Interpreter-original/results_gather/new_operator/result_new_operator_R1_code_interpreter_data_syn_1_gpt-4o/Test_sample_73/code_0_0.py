# Define the operation
def dollar_operation(a, b):
    return 3 * a * b

# Calculate 4 $ 6
result_4_6 = dollar_operation(4, 6)

# Calculate (4 $ 6) $ 4
final_result = dollar_operation(result_4_6, 4)

# Print the final result
print(final_result)