import math

def solve_cube_cutting():
    """
    Calculates the minimum number of cuts to dice a 4x4x4 cube into 1x1x1
    pieces with a knife that has a 2cm cutting depth.
    """
    cube_side = 4
    piece_side = 1
    knife_depth = 2

    # Number of cut planes needed along any single dimension
    cuts_per_dim = cube_side // piece_side - 1

    # --- Dimension 1 ---
    # The first set of cuts is on the original 4cm thick cube.
    thickness_dim1 = cube_side
    # Number of passes needed to cut through the 4cm thickness
    passes_dim1 = math.ceil(thickness_dim1 / knife_depth)
    cuts_dim1 = cuts_per_dim * passes_dim1

    # --- Dimension 2 ---
    # After the first cuts, we have 1cm thick layers. We can stack them.
    # We can create stacks of 2cm height, which matches the knife depth.
    stack_height_dim2 = 2 # e.g., stack two 1cm layers
    passes_dim2 = math.ceil(stack_height_dim2 / knife_depth)
    cuts_dim2 = cuts_per_dim * passes_dim2
    
    # --- Dimension 3 ---
    # After the second cuts, we have 1cm thick rods.
    # We can again create stacks of 2cm height.
    stack_height_dim3 = 2 # e.g., stack two 1cm rods
    passes_dim3 = math.ceil(stack_height_dim3 / knife_depth)
    cuts_dim3 = cuts_per_dim * passes_dim3

    # --- Total and Output ---
    total_cuts = cuts_dim1 + cuts_dim2 + cuts_dim3

    print("Step-by-step calculation for the minimum number of cuts:")
    print("=" * 50)
    
    print(f"To divide a {cube_side}cm length into {piece_side}cm pieces, {cuts_per_dim} cut planes are needed.")
    print("-" * 50)

    # Explanation for Dimension 1
    print("1. Cuts along the first dimension (e.g., Height):")
    print(f"   - The initial cube thickness is {thickness_dim1}cm.")
    print(f"   - With a {knife_depth}cm knife, each cut needs ceil({thickness_dim1}/{knife_depth}) = {passes_dim1} passes.")
    print(f"   - Total cuts for this dimension: {cuts_per_dim} planes * {passes_dim1} passes/plane = {cuts_dim1} cuts.")
    print("-" * 50)
    
    # Explanation for Dimension 2
    print("2. Cuts along the second dimension (e.g., Length):")
    print("   - We now have 1cm thick layers. We can stack them 2 high to make a 2cm stack.")
    print(f"   - The stack height ({stack_height_dim2}cm) allows the {knife_depth}cm knife to cut through in 1 pass.")
    print(f"   - Total cuts for this dimension: {cuts_per_dim} planes * {passes_dim2} pass/plane = {cuts_dim2} cuts.")
    print("-" * 50)
    
    # Explanation for Dimension 3
    print("3. Cuts along the third dimension (e.g., Width):")
    print("   - We now have 1cm thick rods. We can also stack these 2 high.")
    print(f"   - The stack height ({stack_height_dim3}cm) allows the {knife_depth}cm knife to cut through in 1 pass.")
    print(f"   - Total cuts for this dimension: {cuts_per_dim} planes * {passes_dim3} pass/plane = {cuts_dim3} cuts.")
    print("-" * 50)

    # Final Equation and Answer
    print("Total Minimum Cuts Required:")
    print(f"{cuts_dim1} (Dimension 1) + {cuts_dim2} (Dimension 2) + {cuts_dim3} (Dimension 3) = {total_cuts}")
    print("=" * 50)


solve_cube_cutting()