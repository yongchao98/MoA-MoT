def solve_pascal_pyramid():
    """
    Calculates the number of unique values in a 100-layer 3D Square Pascal Pyramid
    based on the provided rules and examples.
    """

    # Step 1: Initialize the pyramid with the layers provided in the problem description.
    # We use these as a starting point due to a slight inconsistency in the generation
    # rule for the very first layers in the user's example.
    layers = [
        [[1]],
        [[1, 1], [1, 1]],
        [[1, 1, 1], [1, 4, 1], [1, 1, 1]],
        [[1, 1, 1, 1], [1, 7, 7, 1], [1, 7, 7, 1], [1, 1, 1, 1]],
        [[1, 1, 1, 1, 1], [1, 10, 16, 10, 1], [1, 16, 28, 16, 1], [1, 10, 16, 10, 1], [1, 1, 1, 1, 1]]
    ]

    # Step 2: Collect all unique numbers from these initial layers using a set.
    unique_numbers = set()
    for layer in layers:
        for row in layer:
            for num in row:
                unique_numbers.add(num)

    # Helper function to get a value from the previous layer.
    # It handles out-of-bounds access by returning 0, as per the pyramid rule.
    def get_value(matrix, r, c):
        rows = len(matrix)
        cols = len(matrix[0])
        if 0 <= r < rows and 0 <= c < cols:
            return matrix[r][c]
        return 0

    # Step 3: Generate the remaining layers from 6 to 100.
    num_layers = 100
    previous_layer = layers[-1]

    for n in range(len(layers) + 1, num_layers + 1):
        # The new layer will have dimensions n x n.
        current_layer = [[0] * n for _ in range(n)]
        
        # Calculate each cell in the new layer.
        for i in range(n):
            for j in range(n):
                # The value is the sum of the 4 cells in the 2x2 block above it.
                val = (get_value(previous_layer, i - 1, j - 1) +
                       get_value(previous_layer, i - 1, j) +
                       get_value(previous_layer, i, j - 1) +
                       get_value(previous_layer, i, j))
                
                current_layer[i][j] = val
                unique_numbers.add(val)
        
        # The newly generated layer becomes the 'previous_layer' for the next iteration.
        previous_layer = current_layer

    # Step 4: Print the final answer, which is the total count of unique numbers.
    print(len(unique_numbers))
    
# Execute the function to find the solution.
solve_pascal_pyramid()