def transform_grid(input_grid):
    # Take the last element and put it at the front
    # Shift everything else right by one
    return [input_grid[-1]] + input_grid[:-1]

# Test examples
examples = [
    ([2,8,0,0,7,0,0,0,0,0,0,0,0,0,0,9,0,2,2,5,0],
     [0,2,8,0,0,7,0,0,0,0,0,0,0,0,0,0,9,0,2,2,5]),
    ([0,3,0,0,0,0,4,2,7,5,3,0,0,0,3,0,0,9,6,0,0],
     [0,0,3,0,0,0,0,4,2,7,5,3,0,0,0,3,0,0,9,6,0]),
    ([5,0,6,4,1,4,6,1,0,0,3,0,8,2,0,0,2,0,6,8,0],
     [0,5,0,6,4,1,4,6,1,0,0,3,0,8,2,0,0,2,0,6,8])
]

# Verify pattern
for ex in examples:
    result = transform_grid(ex[0])
    print(f"Input:   {ex[0]}")
    print(f"Output:  {ex[1]}")
    print(f"Result:  {result}")
    print(f"Matches: {result == ex[1]}\n")

# Apply to test input
test_input = [0,0,0,4,0,9,0,0,3,0,0,3,0,1,0,6,0,0,0,9,0]
result = transform_grid(test_input)
print("Test result:", result)