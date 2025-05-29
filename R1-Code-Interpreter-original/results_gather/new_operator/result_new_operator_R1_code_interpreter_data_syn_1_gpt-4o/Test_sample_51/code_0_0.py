def custom_and(a, b):
    return 2 * a * b

# Calculate (4&4)&4
result = custom_and(custom_and(4, 4), 4)
print(result)