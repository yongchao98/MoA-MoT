# Define the operations
def custom_dollar(a, b):
    return 5 * a + 1 * b

def custom_tilde(a, b):
    return 1 * a * b

def custom_percent(a, b):
    return 2 * a * b

# Calculate the expression (3$4)~5%6
result = custom_percent(custom_tilde(custom_dollar(3, 4), 5), 6)
print(result)