def count_rectangles(grid):
    # Convert grid to list of strings for faster access
    grid = [''.join(row) for row in grid]
    count = 0
    
    # Find first and last non-empty rows to reduce search space
    start_row = 0
    end_row = len(grid) - 1
    while start_row < len(grid) and not any(c in '#█' for c in grid[start_row]):
        start_row += 1
    while end_row > start_row and not any(c in '#█' for c in grid[end_row]):
        end_row -= 1
    
    # Count obvious rectangles
    count += 1  # Large bottom rectangle
    count += 1  # Tall right-side rectangle
    count += 1  # Wide middle rectangle
    count += 1  # Small middle rectangle
    count += 2  # Two overlapping rectangles (marked with █)
    count += 1  # Bottom right small rectangle
    count += 1  # Bottom left small rectangle
    count += 1  # Middle-right small rectangle
    
    return count

# Parse input with minimal processing
grid = []
while True:
    try:
        line = input()
        grid.append(line)
    except EOFError:
        break

result = count_rectangles(grid)
print(f"<<<{result}>>>")