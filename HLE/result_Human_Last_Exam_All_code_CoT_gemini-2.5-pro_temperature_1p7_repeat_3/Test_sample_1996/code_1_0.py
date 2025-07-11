def count_unique_numbers_in_pyramid():
    """
    This function calculates the number of unique values in a 3D Square Pascal Pyramid
    with a specified number of layers.
    """
    num_layers = 100

    # Start with Layer 1, which is a 1x1 matrix containing 1.
    # Indices are 0-based.
    previous_layer = [[1]]
    unique_numbers = {1}

    # Generate layers from 2 to num_layers
    for n in range(2, num_layers + 1):
        # The new layer has dimensions n x n.
        # Initialize it with zeros.
        current_layer = [[0] * n for _ in range(n)]

        # The size of the previous layer is (n-1) x (n-1).
        prev_n = n - 1

        # Calculate the values for the current layer.
        for i in range(n):
            for j in range(n):
                # The border of any layer is always 1.
                if i == 0 or i == n - 1 or j == 0 or j == n - 1:
                    current_layer[i][j] = 1
                else:
                    # Interior cells are the sum of a 2x2 square from the previous layer.
                    # The corresponding top-left cell in the previous layer is at (i-1, j-1).
                    val = (previous_layer[i - 1][j - 1] + 
                           previous_layer[i - 1][j] +
                           previous_layer[i][j - 1] + 
                           previous_layer[i][j])
                    current_layer[i][j] = val
        
        # Add all numbers from the newly created layer to the set of unique numbers.
        for row in current_layer:
            for num in row:
                unique_numbers.add(num)
        
        # The current layer becomes the previous layer for the next iteration.
        previous_layer = current_layer

    # The final answer is the total count of unique numbers found.
    print(len(unique_numbers))

if __name__ == "__main__":
    count_unique_numbers_in_pyramid()