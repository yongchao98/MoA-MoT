# Define inputs and outputs
inputs = [6, 87, 77, 13]
outputs = [8, 170, 150, 22]

# Calculate the differences between input and output
differences = [outputs[i] - inputs[i] for i in range(len(inputs))]

# Check if the differences are related to the input values
ratios = [differences[i] / inputs[i] for i in range(len(inputs))]
print(ratios)