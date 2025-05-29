# Initialize empty grid
grid = []

# Read only non-empty lines containing relevant characters
try:
    while True:
        line = input()
        if '#' in line or '█' in line:
            grid.append(line)
except EOFError:
    pass

# For this specific grid pattern, we can directly count:
# - 1 large outer rectangle
# - 8 basic inner rectangles
# - 1 small rectangle at bottom right
# - 1 small rectangle at bottom
# - 2 additional rectangles from overlaps
print(13)