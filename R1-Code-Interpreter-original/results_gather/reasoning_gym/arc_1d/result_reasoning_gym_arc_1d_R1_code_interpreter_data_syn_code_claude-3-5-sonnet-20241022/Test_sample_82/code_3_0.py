def find_solution(input_str):
    # Convert input to list of integers
    grid = [int(x) for x in input_str.split()]
    n = len(grid)
    output = [0] * n
    
    # Find non-zero sequences
    sequences = []
    start = None
    for i in range(n):
        if grid[i] != 0:
            if start is None:
                start = i
        elif start is not None:
            sequences.append((start, i))
            start = None
    if start is not None:
        sequences.append((start, n))
    
    # Process each sequence
    for start, end in sequences:
        sequence = grid[start:end]
        if start == 0:
            # Rotate right by 4
            rotation = 4
            for i, num in enumerate(sequence):
                new_pos = (i + rotation) % len(sequence)
                output[new_pos] = num
        else:
            # Shift left by 3
            new_start = (start - 3) % n
            for i, num in enumerate(sequence):
                output[(new_start + i) % n] = num
    
    return " ".join(map(str, output))

# Test with the input
test_input = "6 8 4 9 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 9 2"
print(find_solution(test_input))