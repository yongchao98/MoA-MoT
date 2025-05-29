# Initial state
number = 0
light = 'red'  # Initial light state

# Sequence of button presses
sequence = []

# Press A: Add 2, light toggles to green
number += 2
light = 'green'
sequence.append('A')

# Press B: Multiply by 2 (light is green), light toggles to red
number *= 2
light = 'red'
sequence.append('B')

# Press C: Multiply by 2 (light is red), light toggles to green
number *= 2
light = 'green'
sequence.append('C')

# Output the sequence and final number
print("Sequence:", " → ".join(sequence))
print("Final Number:", number)