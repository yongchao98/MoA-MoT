def transform_sequence(grid):
    n = len(grid)
    result = list(grid)
    
    # Find sequences of identical numbers
    i = 0
    while i < n:
        if result[i] != 0:
            start = i
            # Find sequence
            while i < n and result[i] != 0:
                i += 1
            end = i
            
            # Process sequence
            if end - start > 1:  # Need at least 2 numbers
                sequence = result[start:end]
                main_num = None
                diff_num = None
                diff_pos = None
                
                # Find the main number and different number
                for j in range(len(sequence)):
                    if main_num is None:
                        main_num = sequence[j]
                    elif sequence[j] != main_num:
                        diff_num = sequence[j]
                        diff_pos = j
                        break
                    
                if diff_num is None:
                    for j in range(len(sequence)-1, -1, -1):
                        if sequence[j] != main_num:
                            diff_num = sequence[j]
                            diff_pos = j
                            break
                
                if diff_num is not None:
                    # Rearrange based on number comparison
                    if diff_num > main_num:
                        # Move diff_num to front
                        new_sequence = [diff_num] + [main_num] * (len(sequence)-1)
                    else:
                        # Move diff_num to end
                        new_sequence = [main_num] * (len(sequence)-1) + [diff_num]
                    result[start:end] = new_sequence
        i += 1
    return result

# Test input
test_input = [4,3,3,3,3,3,3,3,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
output = transform_sequence(test_input)
print(' '.join(map(str, output)))