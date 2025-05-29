# Initial sequence
sequence = [-2, 4, -10, 24, -58, 140]

# Calculate the differences
def calculate_differences(seq):
    return [seq[i+1] - seq[i] for i in range(len(seq) - 1)]

# Calculate the next term using differences
def next_term(seq):
    differences = [seq]
    while len(differences[-1]) > 1:
        differences.append(calculate_differences(differences[-1]))
    
    # Calculate the next difference
    next_diff = differences[-1][0]
    for i in range(len(differences) - 2, -1, -1):
        next_diff += differences[i][-1]
    
    # Calculate the next term
    return seq[-1] + next_diff

# Find the next term in the sequence
next_value = next_term(sequence)
print(next_value)