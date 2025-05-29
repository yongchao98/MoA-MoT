def transform_grid(input_grid):
    # Convert input string to 2D list
    grid = [[int(x) for x in row.split()] for row in input_grid.strip().split('\n')]
    
    # Flip horizontally (reverse each row)
    result = [row[::-1] for row in grid]
    
    # Convert to string format
    return '\n'.join(' '.join(map(str, row)) for row in result)

# Test input
test_input = """1 3 5
3 5 3
8 8 3"""

print(transform_grid(test_input))