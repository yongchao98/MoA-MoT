# Calculate the number of successful pairs
successful_pairs = [
    (27, 7), (27, 21), (27, 35), (27, 42), (21, 42), (42, 51)
]

# Total number of pairs
total_pairs = 21

# Number of successful pairs
num_successful_pairs = len(successful_pairs)

# Calculate probability
probability = num_successful_pairs / total_pairs
print(probability)