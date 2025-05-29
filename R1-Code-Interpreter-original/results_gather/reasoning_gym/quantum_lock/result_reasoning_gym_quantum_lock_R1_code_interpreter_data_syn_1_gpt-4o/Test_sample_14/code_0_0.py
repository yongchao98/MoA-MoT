# Initial state
number = 0
light = 'red'  # Initial light state
target = 8
sequence = []

# Use Button C to reach the target
while number < target:
    sequence.append('C')
    number += 1
    # Toggle the light
    light = 'green' if light == 'red' else 'red'

# Print the sequence of button presses
print(" → ".join(sequence))