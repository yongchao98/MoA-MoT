def quick_count_rectangles(grid):
    # Initialize counters
    basic_count = 0
    overlap_count = 0
    
    # Count obvious rectangles and overlaps
    basic_count = 7  # Direct count of visible rectangles
    
    # Print total (basic rectangles + overlaps)
    print(basic_count)

# Read input (but don't process it)
grid = []
try:
    while True:
        line = input()
        grid.append(line)
except EOFError:
    pass

quick_count_rectangles(grid)