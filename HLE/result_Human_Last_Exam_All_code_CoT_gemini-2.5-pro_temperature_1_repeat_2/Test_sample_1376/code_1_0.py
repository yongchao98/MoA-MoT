import math

def solve_emperors_engraving():
    """
    Calculates the optimal number of squares (N) and circles (M) to maximize
    the total number of engraved characters (K) from a single sheet of material.
    """
    # 1. Define constants
    mat_dims = [(140, 110), (110, 140)] # Two possible orientations of the material
    sq_s = 10
    circ_s = 40 # Bounding box for the 20cm radius circle
    chars_per_sq = 4
    chars_per_circ = 333 # 999 symbols / 3 symbols per character

    best_N, best_M, max_K = 0, 0, 0

    # Function to calculate max squares for a given number of circles and material size
    def calculate_max_n(M, mat_w, mat_h):
        max_cols = mat_w // circ_s
        max_rows = mat_h // circ_s
        
        if M == 0:
            return (mat_w // sq_s) * (mat_h // sq_s)
        
        if M > max_cols * max_rows:
            return -1 # Impossible to fit

        max_n_for_m = 0
        # Iterate through all possible grid configurations (c x r) that can hold M circles
        for c in range(1, max_cols + 1):
            for r in range(1, max_rows + 1):
                if c * r < M:
                    continue

                grid_w = c * circ_s
                grid_h = r * circ_s

                # Calculate squares in the remaining L-shaped area
                # Place the grid in the top-left corner
                n_bottom = (mat_w // sq_s) * ((mat_h - grid_h) // sq_s)
                n_right = ((mat_w - grid_w) // sq_s) * (grid_h // sq_s)
                
                # Calculate squares in the empty cells within the chosen grid
                empty_cells = c * r - M
                sq_per_cell = (circ_s // sq_s) ** 2
                n_empty = empty_cells * sq_per_cell

                total_n = n_bottom + n_right + n_empty
                if total_n > max_n_for_m:
                    max_n_for_m = total_n
                    
        return max_n_for_m

    # 2. Iterate through all possible numbers of circles (M)
    max_possible_m = (140 // circ_s) * (110 // circ_s)
    for M in range(max_possible_m + 1):
        # For each M, find the best N by checking both material orientations
        n1 = calculate_max_n(M, mat_dims[0][0], mat_dims[0][1])
        n2 = calculate_max_n(M, mat_dims[1][0], mat_dims[1][1])
        current_N = max(n1, n2)

        if current_N < 0:
            continue

        # 3. Calculate total characters (K) and update the best result
        current_K = (current_N * chars_per_sq) + (M * chars_per_circ)

        if current_K > max_K:
            max_K = current_K
            best_N = current_N
            best_M = M
            
    # 4. Print the final answer
    print(f"The optimal combination is N={best_N} squares and M={best_M} circles.")
    print(f"This maximizes the total number of engraved characters at K={max_K}.")
    print(f"The final answer is in the format N:M:K")
    print(f"Result: {best_N}:{best_M}:{max_K}")
    print("\nThe equation for the maximal number of characters is:")
    print(f"{best_N} * {chars_per_sq} + {best_M} * {chars_per_circ} = {max_K}")

solve_emperors_engraving()