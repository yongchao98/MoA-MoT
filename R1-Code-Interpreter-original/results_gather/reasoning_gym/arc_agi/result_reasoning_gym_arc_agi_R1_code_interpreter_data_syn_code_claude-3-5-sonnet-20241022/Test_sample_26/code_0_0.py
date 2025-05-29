def analyze_pattern(input_grid, output_grid):
    rows = len(input_grid)
    cols = len(input_grid[0])
    
    def find_pairs(grid):
        pairs = []
        for i in range(rows-1):
            for j in range(cols-1):
                if grid[i][j] == grid[i][j+1] == grid[i+1][j] == grid[i+1][j+1] == 1:
                    pairs.append((i, j))
        return pairs
    
    def check_transformation(in_pairs, out_pairs):
        for pair in in_pairs:
            i, j = pair
            if output_grid[i][j] == 0:  # If transformed to 0
                # Check if there's another 1-pair above this position
                has_pair_above = False
                for y in range(0, i):
                    if (y, j) in in_pairs:
                        has_pair_above = True
                        break
                print(f"Pair at ({i},{j}): {'changed to 0' if output_grid[i][j] == 0 else 'stayed 1'}, pair above: {has_pair_above}")

    input_pairs = find_pairs(input_grid)
    output_pairs = find_pairs(output_grid)
    check_transformation(input_pairs, output_pairs)

# Example test
example1_input = [
    [8,8,8,8,8,1,1],
    [8,8,8,8,8,1,1],
    [8,8,1,1,8,8,8],
    [8,8,1,1,8,8,8],
    [8,8,8,8,8,8,8]
]

example1_output = [
    [8,8,8,8,8,0,0],
    [8,8,8,8,8,0,0],
    [8,8,1,1,8,8,8],
    [8,8,1,1,8,8,8],
    [8,8,8,8,8,8,8]
]

analyze_pattern(example1_input, example1_output)