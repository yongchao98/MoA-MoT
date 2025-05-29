def analyze_and_transform(input_str):
    # Convert input string to list
    grid = [int(x) for x in input_str.split()]
    
    # Find the sequence of non-zero numbers
    non_zero = []
    for i, num in enumerate(grid):
        if num != 0:
            non_zero.append((num, i))
    
    # Calculate new positions (shift left by 2)
    start_pos = non_zero[0][1]  # First non-zero position
    new_start = max(0, start_pos - 2)  # Shift left by 2
    
    # Create output grid
    output = [0] * len(grid)
    for i, (num, _) in enumerate(non_zero):
        if new_start + i < len(grid):
            output[new_start + i] = num
    
    return ' '.join(map(str, output))

# Test input
test_input = "0 0 0 0 0 0 0 0 0 3 3 3 3 3 3 3 3 0 0 0 0 0 0 0 0 0 0 0 0"
result = analyze_and_transform(test_input)
print(result)