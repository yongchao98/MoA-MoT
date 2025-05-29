# Initial state
number = 0
light = 'red'
sequence = []

# Function to toggle light
def toggle_light(current_light):
    return 'green' if current_light == 'red' else 'red'

# Simulate the process
while number != 16:
    if number < 16:
        if light == 'red':
            # Use A to multiply by 2
            number *= 2
            sequence.append('A')
        else:
            # Use C to add 2
            number += 2
            sequence.append('C')
    # Toggle the light after each button press
    light = toggle_light(light)

# Print the sequence
print(" → ".join(sequence))