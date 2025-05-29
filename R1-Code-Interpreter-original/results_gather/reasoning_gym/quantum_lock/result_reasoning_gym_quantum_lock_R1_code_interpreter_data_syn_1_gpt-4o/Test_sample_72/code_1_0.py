# Initial state
number = 0
light = 'red'  # Initial light state

# Sequence of button presses
sequence = []

# Press button A (Add 2 when red)
if light == 'red':
    number += 2
    light = 'green'
    sequence.append('A')

# Press button C (Multiply by 2 when green)
if light == 'green':
    number *= 2
    light = 'red'
    sequence.append('C')

# Press button B (Multiply by 2 when any)
number *= 2
light = 'green'
sequence.append('B')

# Output the sequence and final number
print("Sequence:", " → ".join(sequence))
print("Final Number:", number)