# Calculate the differences between input and output
inputs = [64, 78, 96, 24]
outputs = [124, 152, 188, 44]

differences = [output - input for input, output in zip(inputs, outputs)]
print(differences)