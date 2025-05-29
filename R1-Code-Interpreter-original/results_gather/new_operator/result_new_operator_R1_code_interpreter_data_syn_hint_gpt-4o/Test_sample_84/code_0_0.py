# Define the operations
def operation_hash(a, b):
    return (a + b) / b

def operation_dollar(a, b):
    return 2 * a + 5 * b

def operation_percent(a, b):
    return (a + b) / b

# Calculate step by step
result_3_hash_4 = operation_hash(3, 4)
result_3_hash_4_dollar_5 = operation_dollar(result_3_hash_4, 5)
final_result = operation_percent(result_3_hash_4_dollar_5, 6)

print(final_result)