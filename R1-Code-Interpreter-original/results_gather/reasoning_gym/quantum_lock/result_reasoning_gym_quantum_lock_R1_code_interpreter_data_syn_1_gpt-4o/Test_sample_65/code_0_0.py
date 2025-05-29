# Initial state
number = 0
sequence = []

# Define button operations
def press_A(n):
    return n - 2

def press_B(n):
    return n + 1

def press_C(n):
    return n * 2

# Simulate the sequence
number = press_B(number)  # B
sequence.append('B')
number = press_C(number)  # C
sequence.append('C')
number = press_C(number)  # C
sequence.append('C')
number = press_C(number)  # C
sequence.append('C')
number = press_C(number)  # C
sequence.append('C')

# Output the result
print("Final number:", number)
print("Sequence:", " → ".join(sequence))