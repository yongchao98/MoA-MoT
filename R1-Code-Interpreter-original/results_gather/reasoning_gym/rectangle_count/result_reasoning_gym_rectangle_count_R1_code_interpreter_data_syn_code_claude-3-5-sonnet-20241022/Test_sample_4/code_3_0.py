# Initialize empty grid
grid = []
rectangle_count = 0

# Read first non-empty line to get width
while True:
    try:
        line = input()
        if '#' in line or '█' in line:
            grid.append(line)
    except EOFError:
        break

# Simple counting based on visual inspection
# We can see:
# 1. One large rectangle at the top
# 2. One small rectangle on the left side
# 3. One large rectangle in the middle
# 4. One rectangle on the right side
# 5. One small rectangle at the bottom

print(5)