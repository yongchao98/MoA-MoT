# Initial position and direction
position = [0, 0]
directions = ["north", "east", "south", "west"]
current_direction_index = 0  # Starting facing "north"

# Instructions
instructions = [
    "turn left", "turn right", "turn around", "turn left",
    "take 4 steps", "turn left", "turn around", "take 3 steps"
]

# Process each instruction
for instruction in instructions:
    if instruction == "turn left":
        current_direction_index = (current_direction_index - 1) % 4
    elif instruction == "turn right":
        current_direction_index = (current_direction_index + 1) % 4
    elif instruction == "turn around":
        current_direction_index = (current_direction_index + 2) % 4
    elif instruction.startswith("take"):
        steps = int(instruction.split()[1])
        if directions[current_direction_index] == "north":
            position[1] += steps
        elif directions[current_direction_index] == "east":
            position[0] += steps
        elif directions[current_direction_index] == "south":
            position[1] -= steps
        elif directions[current_direction_index] == "west":
            position[0] -= steps

# Check if we return to the starting point
result = "Yes" if position == [0, 0] else "No"
print(result)