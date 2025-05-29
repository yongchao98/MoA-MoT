# Calculate the average factor
inputs = [52, 73, 21, 91]
outputs = [100, 142, 38, 178]

factors = [outputs[i] / inputs[i] for i in range(len(inputs))]
average_factor = sum(factors) / len(factors)
print(average_factor)