def parse_grid(grid_str):
    # Parse the grid into a 2D array
    return [[int(x) for x in line.split()] for line in grid_str.strip().split('\n')]

def find_pattern(grid):
    # Count frequencies of numbers (excluding 1s)
    freq = {}
    for row in grid:
        for num in row:
            if num != 1:
                freq[num] = freq.get(num, 0) + 1
    
    # Sort by frequency (descending) and then by value (descending) for ties
    sorted_nums = sorted(freq.items(), key=lambda x: (-x[1], -x[0]))
    
    # Take the top 3 numbers
    return [num for num, _ in sorted_nums[:3]]

# Test input
test_input = """1 1 1 1 1 1 1 1 1 1 1 1 1
1 9 9 1 9 9 1 7 7 1 9 9 1
1 1 1 1 1 1 1 1 1 1 1 1 1
1 9 9 1 7 7 1 9 9 1 2 2 1
1 1 1 1 1 1 1 1 1 1 1 1 1
1 4 4 1 9 9 1 9 9 1 4 4 1
1 1 1 1 1 1 1 1 1 1 1 1 1
1 9 9 1 2 2 1 9 9 1 9 9 1
1 1 1 1 1 1 1 1 1 1 1 1 1
1 9 9 1 9 9 1 9 9 1 7 7 1
1 1 1 1 1 1 1 1 1 1 1 1 1
1 7 7 1 9 9 1 2 2 1 9 9 1
1 1 1 1 1 1 1 1 1 1 1 1 1
1 9 9 1 7 7 1 9 9 1 9 9 1
1 1 1 1 1 1 1 1 1 1 1 1 1"""

grid = parse_grid(test_input)
result = find_pattern(grid)
print('\n'.join(map(str, result)))