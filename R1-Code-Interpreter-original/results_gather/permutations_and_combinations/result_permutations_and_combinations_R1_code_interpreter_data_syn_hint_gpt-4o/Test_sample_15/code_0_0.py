# Define the positions
positions = [None] * 6

# Place A in position 3
positions[2] = 'A'

# Place C in position 2 (adjacent to A and to the left of F)
positions[1] = 'C'

# Place F in position 6 (to the right of C and not in position 5)
positions[5] = 'F'

# Place E and D in positions 4 and 5 (adjacent to each other)
positions[3] = 'E'
positions[4] = 'D'

# Place B in the remaining position 1
positions[0] = 'B'

# Print the final arrangement
print(positions)