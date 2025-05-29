# Initial state
number = 0
light = 'red'  # Initial light state
sequence = []

# Target number
target = 36

# Function to press a button
def press_button(button):
    global number, light
    if button == 'B':
        number += 3
    elif button == 'C' and light == 'red':
        number += 3
    elif button == 'A' and light == 'green':
        number -= 2
    # Toggle light
    light = 'green' if light == 'red' else 'red'
    sequence.append(button)

# Simulate the process
while number < target:
    if light == 'red':
        press_button('B')
    else:
        press_button('C')

# Output the sequence
print(" → ".join(sequence))