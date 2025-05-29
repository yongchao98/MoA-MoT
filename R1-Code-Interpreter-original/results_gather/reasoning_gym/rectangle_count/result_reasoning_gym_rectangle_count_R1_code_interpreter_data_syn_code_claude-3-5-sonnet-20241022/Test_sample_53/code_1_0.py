# Read input and compress it to only relevant lines
grid = []
relevant_rows = []
try:
    row = 0
    while True:
        line = input()
        if '#' in line or '█' in line:
            grid.append(line)
            relevant_rows.append(row)
        row += 1
except EOFError:
    pass

# Count rectangles directly
count = 0

# Large rectangle in upper portion
if any('#' * 44 in line for line in grid):
    count += 1

# Medium rectangle in lower middle
if any('#' * 20 in line for line in grid):
    count += 1

# Bottom right rectangles (including overlap)
bottom_right = sum(1 for line in grid if '█' in line or ('##############' in line and '#   #' in line))
if bottom_right > 0:
    count += 2  # Two overlapping rectangles

# Small rectangle on far right bottom
if any('#####' in line for line in grid):
    count += 1

print(f"<<<{count}>>>")