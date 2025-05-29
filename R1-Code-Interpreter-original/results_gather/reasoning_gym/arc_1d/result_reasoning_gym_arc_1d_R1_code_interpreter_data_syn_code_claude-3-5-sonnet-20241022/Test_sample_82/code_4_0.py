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
            # For sequence at start:
            # Last 4 numbers go to start, remaining numbers go to end
            if len(sequence) > 4:
                last_four = sequence[-4:]
                remaining = sequence[:-4]
                # Place last four at start
                for i, num in enumerate(last_four):
                    output[i] = num
                # Place remaining at end
                for i, num in enumerate(remaining):
                    output[-(len(remaining)-i)] = num
            else:
                # If sequence is 4 or fewer numbers, just keep it at start
                for i, num in enumerate(sequence):
                    output[i] = num
        else:
            # Shift left by 3
            new_start = (start - 3) % n
            for i, num in enumerate(sequence):
                output[(new_start + i) % n] = num
    
    return " ".join(map(str, output))

# Test with the input
test_input = "6 8 4 9 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 9 2"
print(find_solution(test_input))