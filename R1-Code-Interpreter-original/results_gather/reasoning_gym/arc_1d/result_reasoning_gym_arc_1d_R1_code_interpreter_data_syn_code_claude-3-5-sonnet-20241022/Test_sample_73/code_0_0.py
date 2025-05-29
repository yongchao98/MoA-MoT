def analyze_grid(input_grid, output_grid=None):
    # Convert string to list of integers
    if isinstance(input_grid, str):
        input_grid = [int(x) for x in input_grid.split()]
    
    # Find position of 2
    pos_2 = input_grid.index(2)
    
    # Find sequence of non-zero numbers (excluding 2)
    sequence = []
    sequence_start = -1
    sequence_num = -1
    
    for i in range(len(input_grid)):
        if input_grid[i] != 0 and input_grid[i] != 2:
            if sequence_num == -1:
                sequence_num = input_grid[i]
                sequence_start = i
            sequence.append(input_grid[i])
        elif len(sequence) > 0 and input_grid[i] != sequence_num:
            break
    
    # For test case
    if output_grid is None:
        # If 2 is to the right of sequence, move sequence right
        # If 2 is to the left of sequence, move sequence left
        new_grid = [0] * len(input_grid)
        new_grid[pos_2] = 2
        
        sequence_len = len(sequence)
        if pos_2 > sequence_start:
            # Move sequence right by 2
            new_start = sequence_start + 2
        else:
            # Move sequence left by 2
            new_start = sequence_start - 2
            
        for i in range(sequence_len):
            new_grid[new_start + i] = sequence_num
            
        print(' '.join(map(str, new_grid)))

# Test input
test_input = "0 0 0 0 0 0 0 7 7 7 7 7 7 7 7 7 7 7 7 7 0 0 0 0 0 0 2 0"
analyze_grid(test_input)