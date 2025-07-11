def solve_laozi_engraving():
    """
    Calculates the optimal number of circular and squared plates to maximize
    the number of engraved characters, based on a grid-based cutting model.
    """
    # --- Define Constants ---
    # Value (number of characters) for each plate type
    VALUE_CIRCLE = 9999
    VALUE_SQUARE = 360
    
    # Dimensions of the raw material
    MATERIAL_W = 140
    MATERIAL_H = 110
    
    # Dimensions of the bounding box for each plate
    CIRCLE_DIM = 40  # 20cm radius means a 40x40cm bounding box
    SQUARE_DIM = 10  # 10cm side

    # --- Store the best configuration found ---
    # n: number of circles, m: number of squares, k: total characters
    best_config = {'n': 0, 'm': 0, 'k': 0}

    def calculate_best_for_orientation(W, H):
        """
        Tests all possible grid layouts for circles for a given material orientation (W x H)
        and updates the global best_config if a better layout is found.
        """
        nonlocal best_config
        
        # Determine the maximum number of circles that can fit along width and height
        max_circles_w = W // CIRCLE_DIM
        max_circles_h = H // CIRCLE_DIM
        
        # Iterate through all possible grid sizes for the circles
        for num_circles_w in range(max_circles_w + 1):
            for num_circles_h in range(max_circles_h + 1):
                
                # Total circular plates for this grid layout
                n = num_circles_w * num_circles_h
                
                # Calculate the dimensions of the rectangular area consumed by the circle grid
                circles_area_w = num_circles_w * CIRCLE_DIM
                circles_area_h = num_circles_h * CIRCLE_DIM
                
                # By placing the circle grid in a corner, the remaining material
                # forms two rectangles, which can be tiled with squares.
                
                # Rectangle 1
                rem_w1 = W - circles_area_w
                rem_h1 = H
                squares_in_r1 = (rem_w1 // SQUARE_DIM) * (rem_h1 // SQUARE_DIM)
                
                # Rectangle 2
                rem_w2 = circles_area_w
                rem_h2 = H - circles_area_h
                squares_in_r2 = (rem_w2 // SQUARE_DIM) * (rem_h2 // SQUARE_DIM)
                
                # Total squared plates
                m = squares_in_r1 + squares_in_r2
                
                # Calculate the total value (characters) for this configuration
                k = (n * VALUE_CIRCLE) + (m * VALUE_SQUARE)
                
                # If this configuration is the best so far, save it
                if k > best_config['k']:
                    best_config = {'n': n, 'm': m, 'k': k}

    # --- Main Execution ---
    
    # 1. Test with the original material orientation (140x110)
    calculate_best_for_orientation(MATERIAL_W, MATERIAL_H)
    
    # 2. Test with the material rotated (110x140)
    calculate_best_for_orientation(MATERIAL_H, MATERIAL_W)
    
    # 3. Extract final best results
    N = best_config['n']
    M = best_config['m']
    K = best_config['k']
    
    # --- Print the Final Answer ---
    
    # As requested, output each number in the final equation
    value_from_circles = N * VALUE_CIRCLE
    value_from_squares = M * VALUE_SQUARE
    print(f"Final Equation: {N} * {VALUE_CIRCLE} + {M} * {VALUE_SQUARE} = {value_from_circles} + {value_from_squares} = {K}")
    
    # Print the answer in the specified N:M:K format
    print(f"{N}:{M}:{K}")

# Run the solver
solve_laozi_engraving()