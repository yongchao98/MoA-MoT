import collections

def solve_gravitational_wavelet_signatures():
    """
    This script determines the wavelet type and letter index for each of the 9 scalograms.
    The solution is derived by deducing the underlying structure of the 3x3 grid, where:
    1. The wavelet type (M, P, G) is determined by the plot's column.
    2. The letter index (a-i) is determined by the parameters k and n, which are systematically
       mapped to the plot's row and column.
    """

    # Mapping from letter index to the (y, z) pair, as given in the problem.
    # We only need 'y' to identify the index, as each 'y' value from 1-9 is unique.
    index_to_params = {
        'a': (5, 0), 'b': (8, 4), 'c': (1, 1), 'd': (4, 2), 'e': (7, 7),
        'f': (3, 1), 'g': (6, 6), 'h': (9, 1), 'i': (2, 1)
    }

    # Create an inverse map from the 'y' parameter to the letter index for easy lookup.
    y_to_letter = {params[0]: letter for letter, params in index_to_params.items()}
    # y_to_letter will be: {5: 'a', 8: 'b', 1: 'c', 4: 'd', 7: 'e', 3: 'f', 6: 'g', 9: 'h', 2: 'i'}

    # Deduced structure of the 3x3 grid:
    # 1. Columns determine the wavelet type and the parameter 'n' (mass ratio).
    #    - Column 1 (plots 1,4,7): n=1 (equal mass), Wavelet 'M' (Mexican Hat)
    #    - Column 2 (plots 2,5,8): n=2 (unequal mass), Wavelet 'G' (Gabor)
    #    - Column 3 (plots 3,6,9): n=3 (very unequal mass), Wavelet 'P' (Paul)
    col_to_wavelet = {0: 'M', 1: 'G', 2: 'P'}
    
    # 2. Rows and columns together determine the parameter 'k' (mass scale) in a Latin Square pattern.
    #    This is a 2D array representing the grid, where k_map[row][col] gives k.
    k_map = [
        [0, 2, 1],  # Row 0 (Plots 1, 2, 3)
        [1, 0, 2],  # Row 1 (Plots 4, 5, 6)
        [2, 1, 0]   # Row 2 (Plots 7, 8, 9)
    ]

    results = []
    # Iterate through each plot from 1 to 9.
    for plot_num in range(1, 10):
        # Convert 1-based plot number to 0-based row and column indices.
        row = (plot_num - 1) // 3
        col = (plot_num - 1) % 3

        # Determine wavelet type from the column.
        wavelet = col_to_wavelet[col]

        # Determine parameters k and n from the grid position.
        k = k_map[row][col]
        n = col + 1  # Columns 0, 1, 2 correspond to n = 1, 2, 3

        # Calculate the 'y' parameter.
        y = 3 * k + n
        
        # Find the letter index corresponding to this 'y' value.
        letter = y_to_letter[y]

        # Append the two-character code to our results.
        results.append(f"{wavelet}{letter}")
        
    # Format the final output string as requested.
    final_answer = "{" + ", ".join(results) + "}"
    print(final_answer)

solve_gravitational_wavelet_signatures()