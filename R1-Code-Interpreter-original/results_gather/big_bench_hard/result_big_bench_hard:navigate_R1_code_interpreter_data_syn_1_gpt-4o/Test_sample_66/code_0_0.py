# Initial position
x, y = 0, 0

# Step 1: Take 5 steps north
y += 5

# Step 2: Turn around and take 10 steps south
y -= 10

# Step 3: Turn right and take 10 steps west
x -= 10

# Step 4: Turn right and take 10 steps north
y += 10

# Check if back to starting point
is_back_to_start = (x == 0 and y == 0)
print(is_back_to_start)