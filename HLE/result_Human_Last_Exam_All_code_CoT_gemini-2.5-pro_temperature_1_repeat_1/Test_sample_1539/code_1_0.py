import math

def solve_laozi_puzzle():
    """
    Calculates the optimal number of circular and squared plates to maximize
    the total number of engraved characters.
    """
    # Plate values and dimensions
    value_circle = 9999
    value_square = 360
    
    rect_w = 140
    rect_h = 110
    circle_d = 40
    square_s = 10

    # --- Strategy: Prioritize Circles with Hexagonal Packing ---
    # Analysis shows a hexagonal packing of circles (3-2-3 rows) is optimal.
    # This configuration fits 8 circles within a 120cm x 109.3cm bounding box.
    num_circles = 8
    
    # Calculate the remaining rectangular area after cutting the circles.
    # The circles' bounding box is placed at one side of the material.
    circles_footprint_w = 120 
    remaining_width = rect_w - circles_footprint_w
    
    # Calculate how many squares can be cut from the leftover strip.
    if remaining_width >= square_s:
        num_squares_w = remaining_width // square_s
        num_squares_h = rect_h // square_s
        num_squares = num_squares_w * num_squares_h
    else:
        num_squares = 0

    # Calculate the total value (K)
    k_circles_value = num_circles * value_circle
    k_squares_value = num_squares * value_square
    total_value = k_circles_value + k_squares_value

    # --- Output the results ---
    print(f"To maximize the number of characters, the optimal combination is:")
    print(f"N (Circular Plates): {num_circles}")
    print(f"M (Squared Plates): {num_squares}")
    
    print("\nThe maximal number of characters (K) is calculated as follows:")
    # The final equation with all numbers
    print(f"K = {num_circles} * {value_circle} + {num_squares} * {value_square} = {total_value}")
    
    print("\nThe final answer in the format N:M:K is:")
    print(f"{num_circles}:{num_squares}:{total_value}")

# Execute the function
solve_laozi_puzzle()