import numpy as np

def count_unique_numbers_in_pyramid():
    """
    Calculates the number of unique integers in a 100-layer Pascal's square pyramid
    based on the specific rules provided.
    """
    unique_numbers = {1}
    # Layer 1 is a 1x1 matrix with the value 1
    prev_layer = np.array([[1]], dtype=object)

    for n in range(2, 101):  # n represents the size of the layer, from 2x2 to 100x100
        # Layer 2 consists of all 1s.
        if n == 2:
            current_layer = np.ones((2, 2), dtype=object)
            prev_layer = current_layer
            continue

        # Initialize the current layer with its border of 1s.
        current_layer = np.ones((n, n), dtype=object)

        # Calculate the inner elements of the layer.
        # The inner part of an n x n matrix starts at index (1, 1) and ends at (n-2, n-2).
        for i in range(1, n - 1):
            for j in range(1, n - 1):
                # The value is the sum of the four corresponding cells in the previous layer.
                value = prev_layer[i - 1, j - 1] + prev_layer[i - 1, j] + prev_layer[i, j - 1] + prev_layer[i, j]
                current_layer[i, j] = value
                unique_numbers.add(value)

        # The current layer becomes the previous layer for the next iteration.
        prev_layer = current_layer

    return len(unique_numbers)

# Execute the function to find the answer.
final_count = count_unique_numbers_in_pyramid()