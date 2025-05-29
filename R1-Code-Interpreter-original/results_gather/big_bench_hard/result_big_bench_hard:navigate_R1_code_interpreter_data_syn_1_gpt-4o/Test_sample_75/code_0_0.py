# Initial position
position = 0

# Steps
position -= 9  # 9 steps backward
position += 2  # 2 steps forward
position -= 1  # 1 step backward
position += 8  # 8 steps forward

# Check if the final position is the same as the starting point
result = "Yes" if position == 0 else "No"
print(result)