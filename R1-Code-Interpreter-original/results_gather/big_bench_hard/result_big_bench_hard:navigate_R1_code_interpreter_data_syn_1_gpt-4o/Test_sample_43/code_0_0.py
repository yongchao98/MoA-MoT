# Initial position
position = 0

# Steps taken
position += 6  # 6 steps north
position += 8  # 8 steps north
position += 7  # 7 steps north
position -= 6  # 6 steps south (after turning around)

# Check if returned to starting point
returned_to_start = (position == 0)

print("Returned to starting point:", returned_to_start)