# Initial position and direction
position = [0, 0]
direction = [0, 1]  # North

# Steps and turns
steps = [9, 5, 2, 10, 10, 8]
turns = [4]  # Turn right after the fourth step

# Process each step
for i, step in enumerate(steps):
    position[0] += step * direction[0]
    position[1] += step * direction[1]
    if i in turns:
        # Turn right: North -> East, East -> South, South -> West, West -> North
        direction = [direction[1], -direction[0]]

# Check if returned to the starting point
returned_to_start = position == [0, 0]
print("Returned to starting point:", returned_to_start)