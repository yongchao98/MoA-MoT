def find_most_frequent_non_border(grid):
    # Define the blocks based on the pattern observed
    blocks = [
        (1, 1, 3, 3),  # Block 1
        (3, 3, 5, 5),  # Block 2
        (5, 5, 7, 7)   # Block 3
    ]
    
    output = []
    
    for block in blocks:
        r1, c1, r2, c2 = block
        freq = {}
        
        for r in range(r1, r2):
            for c in range(c1, c2):
                num = grid[r][c]
                if num != 1:  # Exclude border numbers
                    if num in freq:
                        freq[num] += 1
                    else:
                        freq[num] = 1
        
        # Find the most frequent number in the block
        most_frequent = max(freq, key=freq.get)
        output.append(most_frequent)
    
    return output

# Test input grid
input_grid = [
    [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
    [1, 9, 9, 1, 9, 9, 1, 7, 7, 1, 9, 9, 1],
    [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
    [1, 9, 9, 1, 7, 7, 1, 9, 9, 1, 2, 2, 1],
    [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
    [1, 4, 4, 1, 9, 9, 1, 9, 9, 1, 4, 4, 1],
    [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
    [1, 9, 9, 1, 2, 2, 1, 9, 9, 1, 9, 9, 1],
    [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
    [1, 9, 9, 1, 9, 9, 1, 9, 9, 1, 7, 7, 1],
    [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
    [1, 7, 7, 1, 9, 9, 1, 2, 2, 1, 9, 9, 1],
    [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
    [1, 9, 9, 1, 7, 7, 1, 9, 9, 1, 9, 9, 1],
    [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
]

output = find_most_frequent_non_border(input_grid)
print(output)