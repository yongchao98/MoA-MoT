# Define the custom modulo operation
def custom_mod(a, b):
    return (a + b) / b

# Calculate (4%3)%8
result = custom_mod(custom_mod(4, 3), 8)
print(result)