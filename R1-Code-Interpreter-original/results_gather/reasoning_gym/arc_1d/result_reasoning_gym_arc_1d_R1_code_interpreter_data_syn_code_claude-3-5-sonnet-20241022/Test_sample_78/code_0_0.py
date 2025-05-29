def find_most_common(sequence):
    if not sequence:
        return 0
    return max(set(sequence), key=sequence.count)

def transform_grid(input_grid):
    # Convert input string to list of integers
    grid = [int(x) for x in input_grid.split()]
    
    # Find continuous sequences and their most common values
    result = grid.copy()
    sequence_start = None
    
    for i in range(len(grid)):
        # Start of a new sequence
        if grid[i] != 0 and sequence_start is None:
            sequence_start = i
        # End of a sequence
        elif (grid[i] == 0 or i == len(grid)-1) and sequence_start is not None:
            end = i if grid[i] == 0 else i + 1
            sequence = grid[sequence_start:end]
            most_common = find_most_common(sequence)
            # Replace all numbers in sequence with most common
            for j in range(sequence_start, end):
                result[j] = most_common
            sequence_start = None
            
    # Print result
    print(' '.join(map(str, result)))

# Test input
test_input = "0 0 0 0 0 0 0 0 1 1 1 9 2 2 4 1 1 1 1 4 1 1 9 1 1 1 1 0 0"
transform_grid(test_input)