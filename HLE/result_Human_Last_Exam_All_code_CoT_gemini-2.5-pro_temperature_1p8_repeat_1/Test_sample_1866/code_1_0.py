import math

def calculate_cuts():
    """
    Calculates the minimum number of cuts to dice a 4x4x4 cube into 1x1x1 cubes
    with a knife that has a cutting depth of 2cm.
    """
    
    # Dimensions of the cube
    initial_dim = 4
    knife_depth = 2

    # --- Step 1: Calculate cuts for the first dimension ---
    # To cut a 4cm dimension into 1cm pieces, 2 cutting operations are needed (4->2, then 2->1).
    # Operation 1 (4cm -> 2cm): Cut the original 4x4x4 cube.
    # The penetration depth is 4cm, so we need ceil(4/2) = 2 passes/cuts.
    cuts_op1_dim1 = math.ceil(initial_dim / knife_depth)
    
    # Operation 2 (2cm -> 1cm): Cut the resulting pieces. They are still 4cm thick
    # in the other dimensions, so penetration depth is still 4cm.
    cuts_op2_dim1 = math.ceil(initial_dim / knife_depth)
    
    cuts_dim1 = cuts_op1_dim1 + cuts_op2_dim1
    print(f"Cuts for the first dimension (e.g., X-axis):")
    print(f"  - The cube is 4cm thick, so cutting it in half requires {cuts_op1_dim1} cuts.")
    print(f"  - The resulting pieces are still in a 4cm thick stack, so the next cut also requires {cuts_op2_dim1} cuts.")
    print(f"  - Total = {cuts_dim1} cuts.\n")

    # After cutting the first dimension, we have slabs that are 1cm thick.
    
    # --- Step 2: Calculate cuts for the second dimension ---
    # We can now arrange the 1cm-thick slabs so the knife penetrates the 1cm side.
    penetration_depth_slim = 1
    
    # Operation 1 (4cm -> 2cm):
    cuts_op1_dim2 = math.ceil(penetration_depth_slim / knife_depth)
    
    # Operation 2 (2cm -> 1cm):
    cuts_op2_dim2 = math.ceil(penetration_depth_slim / knife_depth)
    
    cuts_dim2 = cuts_op1_dim2 + cuts_op2_dim2
    print(f"Cuts for the second dimension (e.g., Y-axis):")
    print(f"  - We now have 1cm thick slabs. We arrange them so the knife penetrates the 1cm side.")
    print(f"  - Cutting the 4cm side in half requires {cuts_op1_dim2} cut.")
    print(f"  - Cutting the resulting 2cm side in half also requires {cuts_op2_dim2} cut.")
    print(f"  - Total = {cuts_dim2} cuts.\n")

    # --- Step 3: Calculate cuts for the third dimension ---
    # The logic is the same for the third dimension, as we now have 1x1 rods.
    cuts_op1_dim3 = math.ceil(penetration_depth_slim / knife_depth)
    cuts_op2_dim3 = math.ceil(penetration_depth_slim / knife_depth)
    cuts_dim3 = cuts_op1_dim3 + cuts_op2_dim3
    print(f"Cuts for the third dimension (e.g., Z-axis):")
    print(f"  - We now have 1cm x 1cm rods. We arrange them so the knife penetrates a 1cm side.")
    print(f"  - This requires {cuts_op1_dim3} + {cuts_op2_dim3} cuts for this dimension.")
    print(f"  - Total = {cuts_dim3} cuts.\n")

    # --- Step 4: Sum the results ---
    total_cuts = cuts_dim1 + cuts_dim2 + cuts_dim3
    print("The final equation for the total minimum number of cuts is:")
    print(f"{cuts_dim1} + {cuts_dim2} + {cuts_dim3} = {total_cuts}")

calculate_cuts()