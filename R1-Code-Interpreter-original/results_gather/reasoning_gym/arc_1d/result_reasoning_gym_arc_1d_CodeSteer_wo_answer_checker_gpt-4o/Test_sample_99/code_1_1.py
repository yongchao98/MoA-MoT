from sympy import symbols, Eq, solve

def detect_transformation(input_grid, output_grid):
    # Create symbolic variables for the input grid
    n = len(input_grid)
    x = symbols(f'x0:{n}')
    
    # Create equations based on the transformation rule
    equations = [Eq(output_grid[i], x[i-1] if i > 0 else 0) for i in range(n)]
    
    # Solve the equations to find the transformation rule
    solution = solve(equations, x)
    
    return solution

def apply_transformation(input_grid, transformation):
    # Apply the transformation to the input grid
    n = len(input_grid)
    output_grid = [0] * n
    for i in range(1, n):
        output_grid[i] = input_grid[i-1]
    return output_grid

# Example input and output grids
input_example = [2, 8, 0, 0, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 9, 0, 2, 2, 5, 0]
output_example = [0, 2, 8, 0, 0, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 9, 0, 2, 2, 5]

# Detect the transformation
transformation = detect_transformation(input_example, output_example)

# Test input grid
test_input = [0, 0, 0, 4, 0, 9, 0, 0, 3, 0, 0, 3, 0, 1, 0, 6, 0, 0, 0, 9, 0]

# Apply the transformation to the test input
test_output = apply_transformation(test_input, transformation)

# Output the result
print("<<<" + " ".join(map(str, test_output)) + ">>>")