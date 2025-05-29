# Define the custom operations
def custom_mod(a, b):
    return (a + b) / b

def custom_tilde(a, b):
    return 1 * a + 1 * b

def custom_hash(a, b):
    return 4 * a + 2 * b

# Calculate the expression (3%4)~5#6
result = custom_hash(custom_tilde(custom_mod(3, 4), 5), 6)
print(result)