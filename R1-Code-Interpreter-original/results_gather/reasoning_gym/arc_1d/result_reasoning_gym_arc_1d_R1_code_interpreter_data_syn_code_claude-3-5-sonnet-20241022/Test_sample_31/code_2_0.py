def transform_grid(grid):
    result = list(grid)
    n = len(grid)
    
    # Find continuous sequences
    i = 0
    while i < n:
        if result[i] != 0:
            # Find the sequence
            start = i
            base_num = result[i]
            sequence_len = 1
            
            # Count how many of the same number
            while i + sequence_len < n and result[i + sequence_len] == base_num:
                sequence_len += 1
                
            # Check for different number at edges
            if i + sequence_len < n and result[i + sequence_len] != 0:
                diff_num = result[i + sequence_len]
                if diff_num > base_num:
                    # Move larger number to front
                    result[start:start+sequence_len+1] = [diff_num] + [base_num] * sequence_len
                else:
                    # Move smaller number to end
                    result[start:start+sequence_len+1] = [base_num] * sequence_len + [diff_num]
            elif i > 0 and result[i-1] != 0:
                diff_num = result[i-1]
                if diff_num > base_num:
                    # Move larger number to front
                    result[i-1:i+sequence_len] = [diff_num] + [base_num] * sequence_len
                else:
                    # Move smaller number to end
                    result[i-1:i+sequence_len] = [base_num] * sequence_len + [diff_num]
            
            i += sequence_len
        else:
            i += 1
            
    return result

# Test input
test_input = [4,3,3,3,3,3,3,3,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
output = transform_grid(test_input)
print(' '.join(map(str, output)))