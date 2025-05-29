# Initial position
position = 0

# Execute the instructions
position += 1  # 1 step right
position += 4  # 4 steps right
position -= 5  # 5 steps left

# Check if the final position is the same as the starting point
result = "Yes" if position == 0 else "No"
print(result)