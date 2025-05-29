def transform_grid(input_grid):
    output = []
    found_first_zero = False
    
    for num in input_grid:
        if num == 0:
            output.append(0)
            found_first_zero = True
        elif num == 3:
            if found_first_zero:
                output.append(1)
            else:
                output.append(2)
    
    # Test with given input
    test_input = [3, 3, 3, 0, 3, 3, 3, 3, 3, 3, 3, 0, 0]
    result = transform_grid(test_input)
    print(' '.join(map(str, result)))

# Run the transformation
test_input = [3, 3, 3, 0, 3, 3, 3, 3, 3, 3, 3, 0, 0]
output = []
found_first_zero = False

for num in test_input:
    if num == 0:
        output.append(0)
        found_first_zero = True
    elif num == 3:
        if found_first_zero:
            output.append(1)
        else:
            output.append(2)

print(' '.join(map(str, output)))