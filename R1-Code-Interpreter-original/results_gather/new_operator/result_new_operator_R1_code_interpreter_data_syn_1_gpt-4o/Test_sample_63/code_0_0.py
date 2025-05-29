def custom_operation(a, b):
    return 2 * a + 1 * b

# Calculate (3#9)#3
result = custom_operation(custom_operation(3, 9), 3)
print(result)