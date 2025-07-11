def solve_cube_cutting_puzzle():
    """
    Calculates and explains the minimum number of cuts to divide a 4x4x4 cube
    into 1x1x1 cubes with a knife depth constraint of 2cm.
    """

    z_axis_cuts = 0
    x_axis_cuts = 0
    y_axis_cuts = 0

    print("--- Puzzle Analysis ---")
    print("Goal: Cut a 4x4x4 cube into 1x1x1 cubes (64 total).")
    print("Constraint: The knife can only cut 2cm deep.")
    print("Strategy: Perform the minimum cuts for each axis (Z, X, Y) by stacking pieces efficiently, while respecting the knife's depth limit.")
    print("-" * 25)

    # Step 1: Z-axis cuts (Horizontal slices)
    print("Step 1: Cuts along the Z-axis")
    print("The cube is 4cm tall, so our first cut must be horizontal to reduce the height.")
    # Cut 1: Halving the main cube
    z_axis_cuts += 1
    print(f"  Cut #{z_axis_cuts}: A horizontal cut at the 2cm mark divides the 4x4x4 cube into two 4x4x2 blocks.")
    # Cut 2: Cutting the two resulting blocks
    z_axis_cuts += 1
    print(f"  Cut #{z_axis_cuts}: Place the two 4x4x2 blocks side-by-side. One more horizontal cut through their middle results in four 4x4x1 slabs.")
    print(f"Total cuts for Z-axis = {z_axis_cuts}")
    print("-" * 25)

    total_cuts_so_far = z_axis_cuts

    # Step 2: X-axis cuts
    print("Step 2: Cuts along the X-axis")
    print("We have four 4x4x1 slabs (1cm high). We can stack two at a time to make a 2cm stack.")
    # Cut 3
    x_axis_cuts += 1
    print(f"  Cut #{total_cuts_so_far + x_axis_cuts}: Make two stacks of two slabs. One vertical cut at the 2cm mark on both stacks creates eight 2x4x1 pieces.")
    # Cut 4
    x_axis_cuts += 1
    print(f"  Cut #{total_cuts_so_far + x_axis_cuts}: Stack the eight 2x4x1 pieces into four 2cm-high stacks. One vertical cut through their middle produces sixteen 1x4x1 rods.")
    print(f"Total cuts for X-axis = {x_axis_cuts}")
    print("-" * 25)

    total_cuts_so_far += x_axis_cuts

    # Step 3: Y-axis cuts
    print("Step 3: Cuts along the Y-axis")
    print("We have sixteen 1x4x1 rods (1cm high). We continue stacking them in pairs.")
    # Cut 5
    y_axis_cuts += 1
    print(f"  Cut #{total_cuts_so_far + y_axis_cuts}: Make eight 2cm-high stacks from the 16 rods. One cut at the 2cm mark on their length creates thirty-two 1x2x1 pieces.")
    # Cut 6
    y_axis_cuts += 1
    print(f"  Cut #{total_cuts_so_far + y_axis_cuts}: Make sixteen 2cm-high stacks. The final cut through their middle yields sixty-four 1x1x1 cubes.")
    print(f"Total cuts for Y-axis = {y_axis_cuts}")
    print("-" * 25)

    # Final Calculation
    total_cuts = z_axis_cuts + x_axis_cuts + y_axis_cuts
    print("Final Result:")
    print("The minimum number of cuts is the sum of cuts for each axis.")
    print(f"{z_axis_cuts} (for Z) + {x_axis_cuts} (for X) + {y_axis_cuts} (for Y) = {total_cuts}")

solve_cube_cutting_puzzle()