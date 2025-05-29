from sympy import symbols, Eq, solve

def transform_grid(input_grid):
    # Define symbolic variables for each position in the grid
    x = symbols('x0:12')
    
    # Create equations based on the input grid
    equations = [Eq(x[i], input_grid[i]) for i in range(len(input_grid))]
    
    # Solve for the positions of zeros
    zero_positions = [i for i, val in enumerate(input_grid) if val == 0]
    non_zero_positions = [i for i, val in enumerate(input_grid) if val != 0]
    
    # Determine the middle position to move zeros
    middle_position = len(input_grid) // 2 - len(zero_positions) // 2
    
    # Create a new grid with zeros in the middle
    output_grid = [0] * len(input_grid)
    for i, pos in enumerate(non_zero_positions):
        if i < middle_position:
            output_grid[i] = input_grid[pos]
        else:
            output_grid[i + len(zero_positions)] = input_grid[pos]
    
    # Place zeros in the middle
    for i in range(middle_position, middle_position + len(zero_positions)):
        output_grid[i] = 0
    
    return output_grid

# Test input
input_grid = [6, 6, 6, 6, 6, 6, 6, 6, 0, 0, 0, 0]
output_grid = transform_grid(input_grid)
print("<<<" + " ".join(map(str, output_grid)) + ">>>")