import math

def solve_cube_cutting():
    """
    Calculates the minimum number of cuts to dice a 4x4x4 cube
    with a knife that has a limited cutting depth of 2cm.
    """
    cube_side = 4
    knife_depth = 2

    # Step 1: Calculate the cuts needed for the three center planes (x=2, y=2, z=2).
    # The cube is 4cm thick, so each center cut requires multiple passes.
    # Passes needed per center plane = ceil(cube_side / knife_depth)
    passes_per_center_plane = math.ceil(cube_side / knife_depth)
    
    num_center_planes = 3
    center_cuts = num_center_planes * passes_per_center_plane
    
    # After the center cuts, we have 2x2x2 blocks.
    
    # Step 2: Calculate cuts for the remaining planes.
    # We now have to cut the 2x2x2 blocks in half along each axis.
    # The thickness is now 2cm, so the 2cm knife needs only one pass.
    # We can stack the pieces to make all cuts for one axis at the same time.
    remaining_cuts_per_axis = 1
    num_axes = 3
    remaining_cuts = num_axes * remaining_cuts_per_axis
    
    # Step 3: Calculate the total cuts.
    total_cuts = center_cuts + remaining_cuts

    print("--- Calculation of Minimum Cuts ---")
    print(f"The cube is {cube_side}x{cube_side}x{cube_side} and the knife can cut {knife_depth}cm deep.")
    print("\nPart 1: The three center cuts")
    print(f"The cube is {cube_side}cm thick, so each center cut requires {passes_per_center_plane} passes.")
    print(f"Number of center cuts = 3 planes * {passes_per_center_plane} passes = {center_cuts} cuts.")
    
    print("\nPart 2: The remaining cuts")
    print("After the first part, we have eight 2x2x2 cubes.")
    print("These are 2cm thick, so they can be cut in a single pass.")
    print("By stacking, we need 1 cut for the X-axis, 1 for the Y-axis, and 1 for the Z-axis.")
    print(f"Number of remaining cuts = 1 (X) + 1 (Y) + 1 (Z) = {remaining_cuts} cuts.")

    print("\n--- Final Answer ---")
    print("The total minimum number of cuts is the sum of the two parts.")
    # The final equation as requested
    print(f"{center_cuts} + {remaining_cuts} = {total_cuts}")

solve_cube_cutting()