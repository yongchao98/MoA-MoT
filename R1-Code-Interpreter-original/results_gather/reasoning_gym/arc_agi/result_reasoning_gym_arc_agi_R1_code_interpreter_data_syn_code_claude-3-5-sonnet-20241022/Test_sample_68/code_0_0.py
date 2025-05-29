def print_transformations(input_grid):
    # Convert string representation to 2D list
    grid = [[int(x) for x in row.split()] for row in input_grid.strip().split('\n')]
    
    # Get dimensions
    rows, cols = len(grid), len(grid[0])
    
    # Function to rotate matrix 90 degrees clockwise
    def rotate_90(matrix):
        return [[matrix[j][i] for j in range(len(matrix)-1,-1,-1)] for i in range(len(matrix[0]))]
    
    # Function to flip matrix horizontally
    def flip_horizontal(matrix):
        return [row[::-1] for row in matrix]
    
    # Function to flip matrix vertically
    def flip_vertical(matrix):
        return matrix[::-1]
    
    # Original
    print("Original:")
    for row in grid:
        print(" ".join(map(str, row)))
    
    # Rotate 90
    rotated = rotate_90(grid)
    print("\nRotated 90:")
    for row in rotated:
        print(" ".join(map(str, row)))
    
    # Flip horizontal
    flipped_h = flip_horizontal(grid)
    print("\nFlipped horizontal:")
    for row in flipped_h:
        print(" ".join(map(str, row)))
    
    # Flip vertical
    flipped_v = flip_vertical(grid)
    print("\nFlipped vertical:")
    for row in flipped_v:
        print(" ".join(map(str, row)))

# Test with Example 2
input_grid = """8 3 3 3
8 3 1 3
8 5 8 5
8 8 5 1"""

print_transformations(input_grid)