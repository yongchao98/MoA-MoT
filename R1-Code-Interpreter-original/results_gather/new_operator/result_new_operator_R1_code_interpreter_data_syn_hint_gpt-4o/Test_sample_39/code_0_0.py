# Define the operation
def custom_operation(a, b):
    return 3 * a * b

# Calculate 7 $ 5
result_7_5 = custom_operation(7, 5)

# Calculate (7 $ 5) $ 7
final_result = custom_operation(result_7_5, 7)

# Print the final result
print(final_result)