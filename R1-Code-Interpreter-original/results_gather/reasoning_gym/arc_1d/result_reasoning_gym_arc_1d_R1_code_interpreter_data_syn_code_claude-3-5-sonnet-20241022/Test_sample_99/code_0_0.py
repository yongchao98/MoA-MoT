def verify_pattern(examples):
    for ex in examples:
        input_grid = ex[0]
        output_grid = ex[1]
        # Perform left rotation
        rotated = input_grid[1:] + [input_grid[0]]
        print(f"Input:  {input_grid}")
        print(f"Output: {output_grid}")
        print(f"Rotated:{rotated}")
        print(f"Matches: {rotated == output_grid}\n")

# Test examples
examples = [
    ([2,8,0,0,7,0,0,0,0,0,0,0,0,0,0,9,0,2,2,5,0],
     [0,2,8,0,0,7,0,0,0,0,0,0,0,0,0,0,9,0,2,2,5]),
    ([0,3,0,0,0,0,4,2,7,5,3,0,0,0,3,0,0,9,6,0,0],
     [0,0,3,0,0,0,0,4,2,7,5,3,0,0,0,3,0,0,9,6,0]),
    ([5,0,6,4,1,4,6,1,0,0,3,0,8,2,0,0,2,0,6,8,0],
     [0,5,0,6,4,1,4,6,1,0,0,3,0,8,2,0,0,2,0,6,8])
]

verify_pattern(examples)

# Now let's apply the rule to the test input
test_input = [0,0,0,4,0,9,0,0,3,0,0,3,0,1,0,6,0,0,0,9,0]
result = test_input[1:] + [test_input[0]]
print("Test result:", result)