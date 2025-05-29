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
            if number == 0:
                number += 2  # Start with an increment to get out of 0
                sequence.append('C')
            else:
                number *= 2
                sequence.append('A')
        else:
            number += 2
            sequence.append('C')
    elif number > 16:
        number -= 2
        sequence.append('B')
    # Toggle the light after each button press
    light = toggle_light(light)

# Print the sequence
print(" → ".join(sequence))