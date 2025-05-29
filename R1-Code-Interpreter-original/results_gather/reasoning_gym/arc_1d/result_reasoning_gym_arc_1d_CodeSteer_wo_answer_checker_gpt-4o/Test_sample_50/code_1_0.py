def transform_grid(input_grid):
    # Step 1: Identify segments
    segments = []
    current_segment = [input_grid[0]]
    
    for num in input_grid[1:]:
        if num == current_segment[-1]:
            current_segment.append(num)
        else:
            segments.append(current_segment)
            current_segment = [num]
    segments.append(current_segment)
    
    # Step 2: Analyze patterns
    non_zero_segments = [seg for seg in segments if seg[0] != 0]
    zero_segments = [seg for seg in segments if seg[0] == 0]
    
    # Step 3: Transform the grid
    # Place the longest non-zero segment at the start
    longest_non_zero_segment = max(non_zero_segments, key=len)
    remaining_non_zero_segments = [seg for seg in non_zero_segments if seg != longest_non_zero_segment]
    
    # Construct the output grid
    output_grid = longest_non_zero_segment + [0] * sum(len(seg) for seg in zero_segments) + sum(remaining_non_zero_segments, [])
    
    return output_grid

# Test input
input_grid = [2, 2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]
output_grid = transform_grid(input_grid)
print("<<<" + " ".join(map(str, output_grid)) + ">>>")