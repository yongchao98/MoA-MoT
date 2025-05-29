# Initial position
position = 0

# Steps
position -= 8  # 8 steps left
position += 6  # 6 steps right
position += 10 # 10 steps right
position -= 2  # 2 steps left
position -= 5  # 5 steps left
position -= 1  # 1 step left

# Check if we return to the starting point
print(position == 0)