def solve_grid(input_str):
    # Convert input to list
    grid = [int(x) for x in input_str.split()]
    
    # Find the majority number
    from collections import Counter
    majority = Counter(grid).most_common(1)[0][0]
    
    # Create output list starting with majority numbers
    output = [majority] * len(grid)
    
    # Collect all positions of non-majority numbers
    non_majority_positions = [i for i, x in enumerate(grid) if x != majority]
    
    # If there are non-majority numbers, shift them all 2 positions left
    if non_majority_positions:
        # Get all non-majority numbers in order
        non_majority_numbers = [grid[i] for i in non_majority_positions]
        
        # Place them starting 2 positions left
        new_start = max(0, non_majority_positions[0] - 2)
        for i, num in enumerate(non_majority_numbers):
            if new_start + i < len(grid):
                output[new_start + i] = num
    
    return ' '.join(map(str, output))

# Test with the given input
test_input = "0 0 0 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 0 0 0 0 0 0 0"
print(solve_grid(test_input))