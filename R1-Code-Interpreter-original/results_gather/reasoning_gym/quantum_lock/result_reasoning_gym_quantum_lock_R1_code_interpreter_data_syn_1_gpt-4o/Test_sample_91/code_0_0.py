# Initial number
number = 0

# Target number
target = 8

# Button B effect
button_B_effect = 2

# Calculate the number of times to press button B
button_B_presses = (target - number) // button_B_effect

# Print the sequence of button presses
sequence = "B " * button_B_presses
print(sequence.strip().replace(" ", " → "))