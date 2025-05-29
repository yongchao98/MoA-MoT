def find_sequence_and_transform(grid):
    n = len(grid)
    # Convert to list for manipulation
    result = list(grid)
    
    # Find sequences of identical numbers
    i = 0
    while i < n:
        # Find start of sequence
        if result[i] != 0:
            start = i
            val = result[i]
            # Find end of sequence
            while i < n and (result[i] == val or (i == start and result[i] != val)):
                i += 1
            end = i
            
            # Check if sequence has different number at start or end
            if start < end-1:  # Sequence must be at least 2 numbers
                if result[start] != val:
                    # Different number at start, move to end
                    diff_num = result[start]
                    result[start:end-1] = [val] * (end-start-1)
                    result[end-1] = diff_num
                elif end < n and result[end-1] != val:
                    # Different number at end, move to start
                    diff_num = result[end-1]
                    result[start+1:end] = [val] * (end-start-1)
                    result[start] = diff_num
        i += 1
    return result

# Test input
test_input = [4,3,3,3,3,3,3,3,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
output = find_sequence_and_transform(test_input)
print(' '.join(map(str, output)))