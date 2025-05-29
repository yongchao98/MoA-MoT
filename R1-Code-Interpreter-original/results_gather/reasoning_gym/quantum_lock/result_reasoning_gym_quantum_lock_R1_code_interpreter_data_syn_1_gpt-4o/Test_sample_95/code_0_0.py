# Initial state
number = 0
light = 'red'  # red or green

# Sequence of button presses
sequence = ['B', 'C', 'A']

# Process each button press
for button in sequence:
    if button == 'A' and light == 'red':
        number *= 2
    elif button == 'B' and light == 'red':
        number += 2
    elif button == 'C' and light == 'green':
        number *= 2
    
    # Toggle the light
    light = 'green' if light == 'red' else 'red'

# Output the final number
print(number)