# Initial state
number = 0
sequence = []

# Step 1: Add 2
number += 2
sequence.append('C')

# Step 2: Multiply by 2
number *= 2
sequence.append('B')

# Step 3: Multiply by 2
number *= 2
sequence.append('B')

# Output the sequence
print(" → ".join(sequence))