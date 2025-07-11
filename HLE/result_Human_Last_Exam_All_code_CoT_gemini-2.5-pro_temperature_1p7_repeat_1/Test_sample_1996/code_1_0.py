import sys

# It is recommended to increase recursion limit for deep calculations, although not strictly necessary for 100 layers.
sys.setrecursionlimit(2000)

def solve():
    """
    Calculates the number of unique values in a 3D Square Pascal Pyramid with 100 layers.
    """
    num_layers = 100

    # A set to store all unique numbers found. Initialize with 1.
    unique_numbers = {1}

    # Layer 1 is just [[1]].
    prev_layer = [[1]]

    # Generate layers from 2 to 100.
    for layer_num in range(2, num_layers + 1):
        size = layer_num
        
        # The new layer has a border of 1s. We can initialize it this way.
        current_layer = [[0] * size for _ in range(size)]

        # Calculate all values for the current layer.
        for r in range(size):
            for c in range(size):
                # The border is always 1.
                if r == 0 or r == size - 1 or c == 0 or c == size - 1:
                    current_layer[r][c] = 1
                else:
                    # Interior cells are the sum of the four cells 'above' in the previous layer.
                    # The four cells are at indices (r-1, c-1), (r-1, c), (r, c-1), and (r, c)
                    # relative to the current cell's position, but in the previous layer's grid.
                    value = (prev_layer[r - 1][c - 1] + 
                             prev_layer[r - 1][c] + 
                             prev_layer[r][c - 1] + 
                             prev_layer[r][c])
                    
                    current_layer[r][c] = value
                    unique_numbers.add(value)

        # The current layer becomes the previous layer for the next iteration.
        prev_layer = current_layer

    print(f"The number of unique numbers in a Pascal Pyramid with {num_layers} layers is:")
    print(len(unique_numbers))

solve()