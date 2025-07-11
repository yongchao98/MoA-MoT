import math

def solve_cutting_problem():
    """
    Calculates the minimum number of cuts to divide a 4x4x4 cube into 1x1x1 cubes
    with a knife that has a cutting depth of 2cm.
    """
    # Problem parameters
    cube_dim = 4
    knife_depth = 2
    num_dims_to_cut = 3

    # We need to slice a 4cm length into four 1cm pieces.
    # The most efficient method is bisection, which involves log2(4) = 2 stages.

    # Stage 1: Halve the 4cm piece into two 2cm pieces.
    thickness_stage1 = 4
    # A single cut through a 4cm piece requires multiple passes.
    # Number of passes = ceiling(thickness / knife_depth)
    cuts_stage1 = math.ceil(thickness_stage1 / knife_depth)

    # Stage 2: Halve the resulting 2cm pieces into 1cm pieces.
    thickness_stage2 = 2
    # After stage 1, we have two 2cm-thick pieces. We stack them to cut simultaneously.
    # The thickness of the pieces being cut is 2cm.
    cuts_stage2 = math.ceil(thickness_stage2 / knife_depth)

    # The total number of cuts for one dimension is the sum of cuts in all stages.
    cuts_per_dim = cuts_stage1 + cuts_stage2

    # This process is repeated for all three dimensions.
    total_cuts = cuts_per_dim * num_dims_to_cut

    print("To solve this, we calculate the cuts for one dimension and then multiply by three.")
    print("A 'cut' is a single pass of the knife.\n")
    print("--- Calculating Cuts for One Dimension (e.g., length) ---")
    print("The 4cm length must be cut into four 1cm slices. This is done in two stages:\n")
    
    print("Stage A: Cut the initial 4cm block into two 2cm blocks.")
    print(f"The block is {thickness_stage1}cm thick, but the knife only cuts {knife_depth}cm deep.")
    print(f"Thus, it requires {cuts_stage1} passes (one from each side) to cut through.")
    
    print("\nStage B: Cut the two 2cm blocks into four 1cm blocks.")
    print("We now have two 2cm-thick blocks. We arrange them side-by-side to cut them at once.")
    print(f"The thickness to cut is {thickness_stage2}cm. A single pass is sufficient.")
    print(f"This takes {cuts_stage2} pass.")

    print(f"\nTotal cuts for one dimension = (Cuts for Stage A) + (Cuts for Stage B) = {cuts_stage1} + {cuts_stage2} = {cuts_per_dim}\n")
    
    print("--- Calculating Total Cuts for the 4x4x4 Cube ---")
    print("The same process is repeated for each of the three dimensions.")
    print(f"Cuts for first dimension = {cuts_per_dim}")
    print(f"Cuts for second dimension = {cuts_per_dim}")
    print(f"Cuts for third dimension = {cuts_per_dim}\n")

    first_dim_cuts_str = str(cuts_per_dim)
    second_dim_cuts_str = str(cuts_per_dim)
    third_dim_cuts_str = str(cuts_per_dim)
    total_cuts_str = str(total_cuts)

    print(f"Total minimum number of cuts = {first_dim_cuts_str} + {second_dim_cuts_str} + {third_dim_cuts_str} = {total_cuts_str}")

solve_cutting_problem()