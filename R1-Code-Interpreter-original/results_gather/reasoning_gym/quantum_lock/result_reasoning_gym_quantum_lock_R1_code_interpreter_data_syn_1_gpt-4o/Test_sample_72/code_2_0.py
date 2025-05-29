# Initial state
number = 0
light = 'red'  # Initial light state

# Sequence of button presses
sequence = ['A', 'C', 'A', 'A']

# Process each button press
for button in sequence:
    if button == 'A' and light == 'red':
        number += 2
    elif button == 'B':
        number *= 2
    elif button == 'C' and light == 'green':
        number *= 2
    
    # Toggle the light
    light = 'green' if light == 'red' else 'red'

# Output the final number
print(number)