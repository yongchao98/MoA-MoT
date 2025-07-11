import numpy as np

def solve_knot_tb():
    """
    Calculates the Thurston-Bennequin number for a given 5x5 grid diagram.
    """
    n = 5
    o_positions = [(1, 1), (2, 2), (3, 3), (4, 4), (5, 5)]
    x_positions = [(1, 4), (2, 5), (3, 1), (4, 2), (5, 3)]

    print(f"The grid size is n = {n}.")
    
    # Create mappings for easy look-up of coordinates
    # For a vertical segment in column i, find the y-coordinates of o and x
    o_y_coords = {pos[0]: pos[1] for pos in o_positions}
    x_y_coords = {pos[0]: pos[1] for pos in x_positions}

    # For a horizontal segment in row j, find the x-coordinates of o and x
    o_x_coords = {pos[1]: pos[0] for pos in o_positions}
    x_x_coords = {pos[1]: pos[0] for pos in x_positions}

    # Determine the orientation sign for each segment based on the path
    # o -> x for horizontal, and x -> o for vertical.
    # Vertical sign s_v[i] = sign(o_y - x_y)
    # Horizontal sign s_h[j] = sign(x_x - o_x)
    s_v = {i: np.sign(o_y_coords[i] - x_y_coords[i]) for i in range(1, n + 1)}
    s_h = {j: np.sign(x_x_coords[j] - o_x_coords[j]) for j in range(1, n + 1)}

    print("\nStep 1: Calculate the writhe (w).")
    print("The writhe is the sum of signs of all crossings.")
    print("A crossing exists at grid point (i, j) if the vertical line in column 'i' intersects the horizontal line in row 'j'.")

    writhe = 0
    crossing_details = []

    for i in range(1, n + 1):  # Column for the vertical segment
        for j in range(1, n + 1):  # Row for the horizontal segment
            
            vy_1, vy_2 = o_y_coords[i], x_y_coords[i]
            hx_1, hx_2 = o_x_coords[j], x_x_coords[j]

            # A crossing exists if the point (i,j) is strictly between the endpoints of both segments.
            is_crossing_i = min(hx_1, hx_2) < i < max(hx_1, hx_2)
            is_crossing_j = min(vy_1, vy_2) < j < max(vy_1, vy_2)

            if is_crossing_i and is_crossing_j:
                sign = s_v[i] * s_h[j]
                writhe += sign
                crossing_details.append(f"  - Crossing at ({i}, {j}): The sign is ({s_v[i]}) * ({s_h[j]}) = {sign}")

    print("\nFound crossings and their signs:")
    if not crossing_details:
        print("  No crossings were found.")
    else:
        for detail in crossing_details:
            print(detail)
            
    print(f"\nThe total writhe is the sum of these signs: w = {writhe}")
    
    print("\nStep 2: Calculate the Thurston-Bennequin number (tb).")
    print("The formula is: tb = w - n")
    
    tb = writhe - n
    
    print("\nFinal Calculation:")
    print(f"tb = {writhe} - {n}")
    print(f"tb = {tb}")

solve_knot_tb()