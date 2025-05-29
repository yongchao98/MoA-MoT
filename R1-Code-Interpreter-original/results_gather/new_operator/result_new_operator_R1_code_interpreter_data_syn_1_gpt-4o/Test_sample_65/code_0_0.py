# Define the operation a#b
def operation(a, b):
    return 3 * a * b

# Calculate 2#9
result_2_9 = operation(2, 9)

# Calculate (2#9)#1
final_result = operation(result_2_9, 1)

# Print the final result
print(final_result)