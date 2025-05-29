def compare_grids(input_grid, output_grid):
    rows = len(input_grid)
    cols = len(input_grid[0])
    changes = []
    
    for i in range(rows):
        for j in range(cols):
            if input_grid[i][j] != output_grid[i][j]:
                changes.append((i, j, input_grid[i][j], output_grid[i][j]))
    
    return changes

def find_pattern(input_str, output_str):
    # Convert string grids to 2D lists
    def str_to_grid(s):
        return [[int(x) for x in row.split()] for row in s.strip().split('\n')]
    
    input_grid = str_to_grid(input_str)
    output_grid = str_to_grid(output_str)
    
    changes = compare_grids(input_grid, output_grid)
    print(f"Found {len(changes)} changes:")
    for change in changes:
        i, j, old, new = change
        print(f"Position ({i}, {j}): {old} -> {new}")
        # Print surrounding values in input grid
        if i > 0 and i < len(input_grid)-1 and j > 0 and j < len(input_grid[0])-1:
            surrounding = [
                input_grid[i-1][j], input_grid[i+1][j],
                input_grid[i][j-1], input_grid[i][j+1]
            ]
            print(f"Surrounding values: {surrounding}")
    
    return changes

# Test with first example (shortened for brevity)
test_input = """5 5 2 3 5 5 2 3 5 5 2 3 2 2 5 3 2 5 5 3 5 2 2
5 5 5 3 5 2 5 3 2 2 5 3 2 2 2 3 5 5 2 3 5 2 2
5 5 5 3 5 5 2 3 5 2 2 3 2 2 5 3 2 2 5 3 2 2 5"""

test_output = """5 8 2 3 5 5 2 3 5 5 2 3 2 2 5 3 2 8 5 3 5 2 2
8 8 5 3 5 2 5 3 2 2 5 3 2 2 2 3 8 8 2 3 5 2 2
5 5 8 3 5 5 2 3 5 2 2 3 2 2 5 3 2 2 8 3 2 2 5"""

find_pattern(test_input, test_output)