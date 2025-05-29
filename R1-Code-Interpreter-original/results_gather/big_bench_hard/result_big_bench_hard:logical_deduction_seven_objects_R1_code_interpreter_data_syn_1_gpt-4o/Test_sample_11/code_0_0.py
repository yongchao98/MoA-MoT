# Initialize the positions with None
positions = [None] * 7

# Place the known positions
positions[6] = 'yellow'  # Yellow is the rightmost
positions[3] = 'black'   # Black is fourth from the left
positions[4] = 'brown'   # Brown is third from the right

# We know the order of green, blue, orange
# and green is to the right of white

# Try placing the remaining books
for i in range(7):
    if positions[i] is None:
        # Try placing white, green, blue, orange
        if 'white' not in positions:
            positions[i] = 'white'
        elif 'green' not in positions:
            positions[i] = 'green'
        elif 'blue' not in positions:
            positions[i] = 'blue'
        elif 'orange' not in positions:
            positions[i] = 'orange'

print(positions)