import numpy as np

def calculate_pascal_pyramid_unique_numbers(num_layers):
    """
    Calculates the number of unique values in a 3D Square Pascal Pyramid.

    Args:
        num_layers: The total number of layers in the pyramid.

    Returns:
        The count of unique numbers across all layers.
    """
    if num_layers == 0:
        return 0
    if num_layers == 1:
        return 1

    # Using a set to store unique numbers efficiently
    unique_numbers = {1}

    # Initialize with the first layer (1-indexed for layer number)
    # Using 0-indexed for matrix access
    previous_layer = np.array([[1]], dtype=object)

    for n in range(2, num_layers + 1):
        # Current layer is n x n (using n for size, which is layer number)
        current_layer = np.ones((n, n), dtype=object)

        # Calculate the interior values
        # The interior exists for layers n > 2
        # Interior indices are from 1 to n-2
        if n > 2:
            for i in range(1, n - 1):
                for j in range(1, n - 1):
                    # Sum of the 2x2 block from the layer above
                    # In 0-based indexing for an n-sized matrix, interior is [1..n-2]
                    # This maps to indices [i-1..i] and [j-1..j] in the (n-1)-sized previous_layer
                    val = (previous_layer[i - 1, j - 1] +
                           previous_layer[i - 1, j] +
                           previous_layer[i, j - 1] +
                           previous_layer[i, j])
                    current_layer[i, j] = val
                    unique_numbers.add(val)
        
        # For n=2, all elements are 1, which are already in unique_numbers
        # so no extra add is needed. current_layer is already all 1s.
        
        # After processing, the current layer becomes the previous for the next iteration
        previous_layer = current_layer

    return len(unique_numbers)

# Calculate for 100 layers
num_unique_elements = calculate_pascal_pyramid_unique_numbers(100)
# print(num_unique_elements)