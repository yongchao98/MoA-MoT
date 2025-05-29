# Initial state
number = 0
sequence = []

# Step 1: Add 2
number += 2
sequence.append('B')

# Step 2: Multiply by 2
number *= 2
sequence.append('C')

# Step 3: Multiply by 2
number *= 2
sequence.append('C')

# Output the sequence and final number
print("Sequence:", " → ".join(sequence))
print("Final Number:", number)