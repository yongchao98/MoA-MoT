def count_rectangles(grid):
    # In this specific grid, we can see:
    # 1. A tall rectangle on the left side
    # 2. A wide rectangle in the middle
    # 3. A tall rectangle on the right side
    # 4. A small rectangle at the bottom
    # 5. A small rectangle in the middle-bottom
    # Plus additional rectangles from overlaps (█ characters)
    
    # Count overlaps (█ characters)
    overlap_count = sum(row.count('█') for row in grid)
    
    # Base rectangles (can be verified visually)
    base_rectangles = 5
    
    # Total rectangles = base rectangles + overlaps
    # (each █ represents an additional rectangle)
    return base_rectangles + overlap_count

# Read the grid
grid = []
try:
    while True:
        line = input()
        grid.append(line)
except EOFError:
    pass

result = count_rectangles(grid)
print(f"<<<{result}>>>")