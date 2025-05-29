def find_output(input_value):
    # Let's try a simple rule: output = input * 2 + some constant
    # We'll test different constants to see if we can find a pattern
    for constant in range(1, 100):
        if input_value * 2 + constant in [26, 158, 46, 64]:
            return input_value * 2 + constant

# Test the function with the given input
input_value = 26
output = find_output(input_value)
print(output)