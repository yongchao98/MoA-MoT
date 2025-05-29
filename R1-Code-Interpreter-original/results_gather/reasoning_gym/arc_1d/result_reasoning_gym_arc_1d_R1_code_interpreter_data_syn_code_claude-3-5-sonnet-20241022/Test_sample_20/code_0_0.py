def process_grid(input_grid):
    # Count occurrences
    count_2 = input_grid.count(2)
    count_0 = input_grid.count(0)
    count_1 = input_grid.count(1)
    
    # Create output according to rule
    output = [2] * count_2 + [0] * count_0 + [1] * count_1
    
    # Convert to string format
    return ' '.join(map(str, output))

# Test input
test_input = [0, 1, 2, 1, 0, 0, 1, 0, 1, 0, 2]
print(process_grid(test_input))