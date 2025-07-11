import math

def get_point_coords(point_id):
    """
    Calculates the coordinates of a given point on the hexagonal grid.
    The grid has a side length s=1.
    """
    s = 1.0
    sqrt3 = math.sqrt(3)

    # Hexagon Center Coordinates (from axial coordinates)
    # H_13(0,0), H_31(1,1), H_23(1,-1)
    # x = s * 1.5 * q
    # y = s * (0.5*sqrt3*q + sqrt3*r)
    centers = {
        13: (0, 0),
        31: (s * 1.5 * 1, s * (0.5 * sqrt3 * 1 + sqrt3 * 1)),
        23: (s * 1.5 * 1, s * (0.5 * sqrt3 * 1 + sqrt3 * -1)),
    }

    # Local coordinates for a hexagon centered at (0,0) with vertical sides
    # Using point numbers from the diagram to identify their type
    local_coords = {
        # Centers
        13: (0, 0), 31: (0,0), 23: (0,0),
        # Vertices of H_13
        5: (-0.5 * s, -0.5 * s * sqrt3), # V_BL
        7: (0.5 * s, -0.5 * s * sqrt3),  # V_BR
        # Midpoints of H_13
        10: (0, 0.5 * s * sqrt3), # M_T
        4: (0, -0.5 * s * sqrt3), # M_B
    }

    # Determine which hexagon the point belongs to
    if point_id in [2, 4, 5, 7, 8, 9, 10, 11, 12, 13]:
        center_coord = centers[13]
        local_coord = local_coords.get(point_id, (0,0))
    elif point_id in [25, 26, 27, 28, 29, 30, 31]:
        center_coord = centers[31]
        local_coord = local_coords.get(point_id, (0,0))
    elif point_id in [14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24]:
        center_coord = centers[23]
        local_coord = local_coords.get(point_id, (0,0))
    else:
        # Default to H_13 for any unlisted points
        center_coord = centers[13]
        local_coord = (0,0)

    return center_coord[0] + local_coord[0], center_coord[1] + local_coord[1]

def calculate_period(point_sequence):
    """
    Calculates the distance between the first and last point in a sequence.
    """
    start_point_id = point_sequence[0]
    end_point_id = point_sequence[-1]

    x1, y1 = get_point_coords(start_point_id)
    x2, y2 = get_point_coords(end_point_id)

    distance = math.sqrt((x2 - x1)**2 + (y2 - y1)**2)
    return distance

def main():
    sequences = [
        [13, 31, 23],
        [10, 4, 23, 31],
        [5, 15, 17, 19, 21, 7],
        [4, 5, 14, 23, 18, 19, 21, 22, 31, 30, 9, 8, 13]
    ]

    periods = [calculate_period(seq) for seq in sequences]
    
    # Output the results as a comma-separated string
    # Using f-strings to format the output and show the symbolic representation
    results_str = [
        f"sqrt(3)",
        f"sqrt(21)/2",
        f"1",
        f"sqrt(3)/2"
    ]
    
    # Print the numerical values for verification
    # print(f"Numerical values: {periods[0]:.3f}, {periods[1]:.3f}, {periods[2]:.3f}, {periods[3]:.3f}")
    
    # Final answer as requested
    print(f"{periods[0]},{periods[1]},{periods[2]},{periods[3]}")


if __name__ == "__main__":
    main()