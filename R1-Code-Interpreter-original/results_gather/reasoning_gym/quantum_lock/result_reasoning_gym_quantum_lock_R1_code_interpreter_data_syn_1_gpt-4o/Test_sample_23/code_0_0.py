# Initial number and target
start = 0
target = 9

# Sequence of button presses
sequence = ['B', 'B', 'B', 'B', 'A']

# Apply the sequence to the start number
current = start
for button in sequence:
    if button == 'A':
        current += 1
    elif button == 'B':
        current += 2
    elif button == 'C':
        current -= 1

# Print the final result
print(current)