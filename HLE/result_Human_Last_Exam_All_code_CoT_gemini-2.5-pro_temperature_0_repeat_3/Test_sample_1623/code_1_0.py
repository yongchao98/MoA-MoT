import numpy as np

def solve_knot_tb():
    """
    Calculates the maximal Thurston-Bennequin number for a knot from a grid diagram.
    """
    n = 5
    o_coords = [(1, 1), (2, 2), (3, 3), (4, 4), (5, 5)]
    x_coords = [(1, 4), (2, 5), (3, 1), (4, 2), (5, 3)]

    # Create dictionaries for easy lookup
    # For a given column i, find the row of O and X
    o_rows = {c: r for c, r in o_coords}
    x_rows = {c: r for c, r in x_coords}
    # For a given row j, find the column of O and X
    o_cols = {r: c for c, r in o_coords}
    x_cols = {r: c for c, r in x_coords}

    writhe = 0
    crossings = []

    # Iterate over all possible crossing points (i, j)
    for i in range(1, n + 1):
        for j in range(1, n + 1):
            # A crossing cannot be at the same location as an O or X
            if (i, j) in o_coords or (i, j) in x_coords:
                continue

            # Get O and X positions for the current column and row
            o_r_i, x_r_i = o_rows[i], x_rows[i]
            o_c_j, x_c_j = o_cols[j], x_cols[j]

            # Check for crossing condition
            # Vertical segment in column i must span row j
            cond1 = (min(o_r_i, x_r_i) < j < max(o_r_i, x_r_i))
            # Horizontal segment in row j must span column i
            cond2 = (min(o_c_j, x_c_j) < i < max(o_c_j, x_c_j))

            if cond1 and cond2:
                # Calculate the sign of the crossing
                vert_dir = x_r_i - o_r_i
                horiz_dir = x_c_j - o_c_j
                sign = np.sign(vert_dir * horiz_dir)
                
                crossings.append(((i, j), sign))
                writhe += sign

    print("Calculating the writhe by finding all crossings:")
    if not crossings:
        print("No crossings found.")
    else:
        for cross in crossings:
            print(f"Crossing at {cross[0]} has sign {cross[1]}.")
    
    print("\nThe total writhe is the sum of the signs of all crossings.")
    
    # Constructing the equation string for writhe
    writhe_eq_parts = [str(c[1]) for c in crossings]
    if not writhe_eq_parts:
        writhe_eq_str = "0"
    else:
        writhe_eq_str = " + ".join(f"({p})" for p in writhe_eq_parts)

    print(f"Total writhe = {writhe_eq_str} = {writhe}")

    tb = writhe - n
    
    print(f"\nThe grid number n is {n}.")
    print("The maximal Thurston-Bennequin number is calculated as writhe - n.")
    print(f"The maximal Thurston-Bennequin number is {writhe} - {n} = {tb}")

    # The final answer in the required format
    print(f"\n<<<{tb}>>>")

solve_knot_tb()