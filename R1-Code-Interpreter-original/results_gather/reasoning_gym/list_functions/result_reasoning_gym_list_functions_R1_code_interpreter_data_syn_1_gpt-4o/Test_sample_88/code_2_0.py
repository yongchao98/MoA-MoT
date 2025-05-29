def transform(input_value):
    # Let's try a different transformation
    # We will test a multiplication and addition
    factor = 2  # Hypothetical multiplication factor
    addition = 0  # Hypothetical addition
    return input_value * factor + addition

# Test the function with the given examples
examples = [25, 29, 17, 58]
outputs = [46, 54, 30, 112]

# Check if the transformation works for all examples
for i, example in enumerate(examples):
    assert transform(example) == outputs[i], f"Failed for input {example}"

# Apply the transformation to the new input
new_input = 1
new_output = transform(new_input)
print(new_output)