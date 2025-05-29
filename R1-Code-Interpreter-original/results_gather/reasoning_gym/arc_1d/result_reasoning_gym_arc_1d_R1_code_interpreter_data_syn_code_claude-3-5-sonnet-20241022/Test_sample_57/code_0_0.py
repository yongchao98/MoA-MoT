def transform_grid(input_str):
    # Convert input string to list of integers
    grid = [int(x) for x in input_str.split()]
    
    # Find the sequence of non-zero numbers
    non_zero = []
    start_idx = -1
    for i, num in enumerate(grid):
        if num != 0 and start_idx == -1:
            start_idx = i
        if num != 0:
            non_zero.append(num)
    
    # Calculate new starting position
    seq_length = len(non_zero)
    new_start = 8 - seq_length//2
    
    # Create output grid
    output = [0] * 29
    for i in range(seq_length):
        output[new_start + i] = non_zero[i]
    
    return ' '.join(map(str, output))

# Test input
test = "0 0 0 0 0 0 0 0 0 0 0 0 0 6 6 6 6 6 6 6 0 0 0 0 0 0 0 0 0"
print(transform_grid(test))