def solve_braid_index():
    """
    Calculates the braid index of a knot from its grid diagram representation.
    """
    # Grid number
    n = 7

    # Coordinates of the 'o' and 'x' markers
    # (column, row) format, where (1,1) is bottom-left
    o_coords = [(1, 1), (2, 7), (3, 4), (4, 5), (5, 3), (6, 6), (7, 2)]
    x_coords = [(1, 2), (2, 6), (3, 3), (4, 1), (5, 7), (6, 5), (7, 4)]

    # Create dictionaries for quick lookups of 'x' positions
    # x_by_row[j] gives the column i of the 'x' in row j
    x_by_row = {j: i for i, j in x_coords}
    # x_by_col[i] gives the row j of the 'x' in column i
    x_by_col = {i: j for i, j in x_coords}

    # Initialize a counter for North-East (NE) type markers
    ne_count = 0

    print("Analyzing 'o' markers based on their corresponding 'x' positions:")
    print("-" * 60)
    # Loop through each 'o' to classify it
    for i, j in o_coords:
        # Find the 'x' in the same row
        i_x = x_by_row[j]
        # Find the 'x' in the same column
        j_x = x_by_col[i]

        # Classify the 'o' marker
        marker_type = ""
        if i_x > i and j_x > j:
            marker_type = "NE (North-East)"
            ne_count += 1
        elif i_x < i and j_x > j:
            marker_type = "NW (North-West)"
        elif i_x > i and j_x < j:
            marker_type = "SE (South-East)"
        elif i_x < i and j_x < j:
            marker_type = "SW (South-West)"
        
        print(f"o at ({i},{j}): row x at ({i_x},{j}), col x at ({i},{j_x}). Type: {marker_type}")

    print("-" * 60)
    
    # Calculate the braid index using the formula: n - (number of NE markers)
    braid_index = n - ne_count

    print(f"\nGrid number (n): {n}")
    print(f"Number of NE-type 'o' markers: {ne_count}")
    print("\nThe braid index is calculated as: n - (number of NE markers)")
    print(f"Braid Index = {n} - {ne_count} = {braid_index}")

    return braid_index

if __name__ == "__main__":
    result = solve_braid_index()
    print(f"\n<<< {result} >>>")
