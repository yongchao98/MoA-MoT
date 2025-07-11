import collections

def count_unique_numbers_in_pyramid():
    """
    This function calculates the number of unique integers in a 3D Square Pascal Pyramid
    of 100 layers.
    """
    num_layers = 100

    # Layer 1 is a special case, it just contains 1.
    if num_layers == 0:
        print(0)
        return
        
    previous_layer = [[1]]
    unique_numbers = {1}

    # Generate layers from 2 to num_layers
    for n in range(2, num_layers + 1):
        # The current layer has dimensions n x n
        current_layer = [[0] * n for _ in range(n)]
        
        # We can optimize by calculating only the upper triangle (where col >= row)
        # and copying the values to the lower triangle, because the pyramid is symmetric.
        for r in range(n):
            for c in range(r, n):
                # The border cells of any layer are always 1
                if r == 0 or r == n - 1 or c == 0 or c == n - 1:
                    value = 1
                else:
                    # An inner cell's value is the sum of the four cells above it.
                    # Note: The coordinates for the previous layer are different because it's a smaller grid.
                    # The 2x2 block from the previous layer corresponds to:
                    # [r-1, c-1], [r-1, c]
                    # [r, c-1],   [r, c]
                    val1 = previous_layer[r - 1][c - 1]
                    val2 = previous_layer[r - 1][c]
                    val3 = previous_layer[r][c - 1]
                    val4 = previous_layer[r][c]
                    value = val1 + val2 + val3 + val4

                current_layer[r][c] = value
                current_layer[c][r] = value # Apply symmetry
                unique_numbers.add(value)
        
        # The newly generated layer becomes the 'previous_layer' for the next iteration
        previous_layer = current_layer

    print(f"The total number of unique numbers in a {num_layers}-layer pyramid is: {len(unique_numbers)}")

count_unique_numbers_in_pyramid()
<<<2553>>>