def analyze_and_transform(input_grid):
    # Count initial non-zero numbers
    initial_count = 0
    num = None
    for x in input_grid:
        if x == 0:
            break
        if num is None:
            num = x
        initial_count += 1
    
    # Count total non-zeros and zeros
    total_num = sum(1 for x in input_grid if x == num)
    total_zeros = input_grid.count(0)
    
    # Create output grid
    output = []
    # Add initial number + 2 more
    output.extend([num] * (initial_count + 2))
    # Add all zeros
    output.extend([0] * total_zeros)
    # Fill remaining with the number
    remaining = 24 - len(output)
    output.extend([num] * remaining)
    
    return output

# Test input
test = [4, 4, 4, 4, 4, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4]
result = analyze_and_transform(test)
print(' '.join(map(str, result)))