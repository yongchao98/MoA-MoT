import math

def calculate_cuts():
    """
    Calculates the minimum number of cuts to dice a 4x4x4 cube into 1x1x1 cubes
    with a knife that can only cut 2cm deep.
    """
    
    knife_depth = 2
    
    print("To cut a 4x4x4 cube into 1x1x1 cubes, we need 3 cuts along each of the x, y, and z axes.")
    print("The cost of a cut depends on the height of the stacked pieces: cost = ceil(Stack_Height / 2).\n")

    # --- Stage 1: Cut the 4x4x4 cube into eight 2x2x2 cubes ---
    print("--- Stage 1: Cutting the 4x4x4 cube into 2x2x2 cubes ---")
    print("This requires making the three central cuts at the 2cm mark of each axis.")
    
    # For each central cut, we start with a 4x4x4 block (or a stack that is 4cm high).
    stage1_stack_height = 4
    cuts_per_plane_s1 = math.ceil(stage1_stack_height / knife_depth)
    print(f"To make a central cut, the pieces must be stacked, resulting in a height of {stage1_stack_height}cm.")
    print(f"Number of cuts for one central plane = ceil({stage1_stack_height} / {knife_depth}) = {cuts_per_plane_s1} cuts.")
    
    num_central_planes = 3
    total_cuts_s1 = num_central_planes * cuts_per_plane_s1
    print(f"There are {num_central_planes} central planes (x=2, y=2, z=2).")
    print(f"Total cuts for Stage 1 = {num_central_planes} planes * {cuts_per_plane_s1} cuts/plane = {total_cuts_s1} cuts.\n")
    
    # --- Stage 2: Cut the eight 2x2x2 cubes into 1x1x1 cubes ---
    print("--- Stage 2: Cutting the eight 2x2x2 cubes into 1x1x1 cubes ---")
    print("This requires cutting all resulting 2x2x2 pieces in half along each of their three axes.")
    
    # The pieces are now 2x2x2. We can arrange them so the stack height is 2cm.
    stage2_stack_height = 2
    cuts_per_plane_s2 = math.ceil(stage2_stack_height / knife_depth)
    print(f"For this stage, we can arrange the 2x2x2 pieces so the stack height is {stage2_stack_height}cm.")
    print(f"Number of cuts for one plane = ceil({stage2_stack_height} / {knife_depth}) = {cuts_per_plane_s2} cut.")

    num_planes_s2 = 3
    total_cuts_s2 = num_planes_s2 * cuts_per_plane_s2
    print(f"This must be done for all {num_planes_s2} axes (x, y, z).")
    print(f"Total cuts for Stage 2 = {num_planes_s2} planes * {cuts_per_plane_s2} cut/plane = {total_cuts_s2} cuts.\n")

    # --- Final Calculation ---
    print("--- Total Minimum Cuts ---")
    total_cuts = total_cuts_s1 + total_cuts_s2
    print(f"The total minimum number of cuts is the sum of the cuts from both stages.")
    print(f"{total_cuts_s1} + {total_cuts_s2} = {total_cuts}")

calculate_cuts()