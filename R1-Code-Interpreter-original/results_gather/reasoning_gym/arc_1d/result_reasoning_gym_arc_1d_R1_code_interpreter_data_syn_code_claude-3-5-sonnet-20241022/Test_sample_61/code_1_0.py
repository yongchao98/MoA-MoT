def apply_rule(input_grid):
    grid = [int(x) for x in input_grid.split()]
    result = grid.copy()
    
    i = 0
    while i < len(grid):
        if grid[i] == 2:
            start = i
            # Find end of sequence
            while i < len(grid) and grid[i] == 2:
                i += 1
            end = i
            
            seq_length = end - start
            if seq_length <= 3:
                # Convert first number to 8
                result[start] = 8
            else:
                # Convert first seven numbers to 8 (or all but one if shorter)
                for j in range(start, min(start + 7, end)):
                    result[j] = 8
        else:
            i += 1
    
    return ' '.join(map(str, result))

test_input = "2 2 2 2 2 2 2 2 0 2 2 2 2 2 2 2 2 2 0 0 0 0 0 0 0 0 0 0"
print(apply_rule(test_input))