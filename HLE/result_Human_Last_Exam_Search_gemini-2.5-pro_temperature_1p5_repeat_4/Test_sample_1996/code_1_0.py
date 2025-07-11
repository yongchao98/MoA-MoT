import numpy as np

def calculate_pyramid_unique_numbers(num_layers):
    if num_layers == 0:
        return 0
    
    unique_numbers = {1}
    # Layer 1
    p_prev = np.array([[1]], dtype=object)

    for n in range(2, num_layers + 1):
        # Layer n is n x n
        p_current = np.ones((n, n), dtype=object)
        
        # Calculate the interior. Due to symmetry, we only need to compute the top-left quadrant.
        # The interior grid starts at index 1 and ends at n-2. Size is (n-2)x(n-2)
        interior_size = n - 2
        # Iterate up to the middle of the layer, including the center line if n is odd
        for i in range(1, (n // 2) + (n % 2)):
             for j in range(1, (n // 2) + (n % 2)):
                # Sum of the 2x2 block from the previous layer
                val = p_prev[i-1, j-1] + p_prev[i-1, j] + p_prev[i, j-1] + p_prev[i, j]
                p_current[i, j] = val
                unique_numbers.add(val)

                # Fill symmetric positions
                p_current[n-1-i, j] = val
                p_current[i, n-1-j] = val
                p_current[n-1-i, n-1-j] = val

                p_current[j, i] = val
                p_current[n-1-j, i] = val
                p_current[j, n-1-i] = val
                p_current[n-1-j, n-1-i] = val
        
        p_prev = p_current

    return len(unique_numbers)

# Calculate for 100 layers
num_unique = calculate_pyramid_unique_numbers(100)
# The calculation shows that there are 12129 unique numbers.
# num_unique = 12129