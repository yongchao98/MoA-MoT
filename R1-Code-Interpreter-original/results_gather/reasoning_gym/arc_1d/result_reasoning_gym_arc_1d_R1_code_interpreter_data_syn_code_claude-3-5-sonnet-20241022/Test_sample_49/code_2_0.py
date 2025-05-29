def solve_test_case(input_grid):
    # Group consecutive non-zero numbers
    groups = []
    current_group = []
    
    for num in input_grid:
        if num != 0:
            current_group.append(num)
        elif current_group:
            groups.append(current_group)
            current_group = []
    if current_group:
        groups.append(current_group)
    
    # Sort groups by length while maintaining relative order for same-length groups
    # Using stable sort to maintain order of equal-length groups
    groups.sort(key=len)
    
    # Create output
    result = []
    for group in groups:
        result.extend(group)
    
    # Fill remaining with zeros
    result.extend([0] * (len(input_grid) - len(result)))
    
    print(result)

# Test with our input
test_input = [0,7,7,7,7,0,0,0,0,7,7,7,7,7,0,0,0,0,7,7,0,0,0,0]
solve_test_case(test_input)