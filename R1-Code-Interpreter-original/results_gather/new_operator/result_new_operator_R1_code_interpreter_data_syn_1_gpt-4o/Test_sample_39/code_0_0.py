# Define the operation
def dollar_operation(a, b):
    return 3 * a * b

# Calculate 7 $ 5
result_7_5 = dollar_operation(7, 5)

# Calculate (7 $ 5) $ 7
final_result = dollar_operation(result_7_5, 7)

# Print the final result
print(final_result)