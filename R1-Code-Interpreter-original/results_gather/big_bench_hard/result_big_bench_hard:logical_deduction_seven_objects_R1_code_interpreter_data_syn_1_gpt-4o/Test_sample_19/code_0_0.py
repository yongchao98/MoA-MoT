# Initialize the positions with known values
positions = [None] * 7
positions[6] = 'white'  # The white book is the rightmost
positions[4] = 'gray'   # The gray book is third from the right
positions[2] = 'orange' # The orange book is second from the left
positions[3] = 'black'  # The black book is third from the left

# Determine the positions of the remaining books
# The green book is to the left of the gray book
# The yellow book is to the left of the black book

# Possible positions for green and yellow are 0 and 1
# Since orange is at 2 and black is at 3

# Try placing green and yellow
if positions[0] is None:
    positions[0] = 'green'
    positions[1] = 'yellow'
else:
    positions[0] = 'yellow'
    positions[1] = 'green'

# The remaining position is for the purple book
for i in range(7):
    if positions[i] is None:
        positions[i] = 'purple'

# Output the book in the third position from the left
print(positions[2])