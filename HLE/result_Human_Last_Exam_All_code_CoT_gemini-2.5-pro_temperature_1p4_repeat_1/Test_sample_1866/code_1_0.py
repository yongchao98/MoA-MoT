import math

def solve_cube_cutting():
    """
    Calculates the minimum number of cuts to get 1x1x1 cubes from a 4x4x4 cube
    with a knife that has a 2cm cutting depth.
    """
    
    # --- Stage 1: Cut the 4x4x4 cube into eight 2x2x2 cubes ---
    # This requires making one central cut along each of the three axes.
    # For each cut, we arrange the pieces to form a 4cm high block.
    # The number of passes for one such cut is ceil(4cm / 2cm) = 2.
    
    passes_per_4cm_cut = math.ceil(4 / 2)
    
    # Cut 1 (z-axis): Cut the 4x4x4 cube at the 2cm mark.
    stage1_z_cuts = passes_per_4cm_cut
    
    # Cut 2 (y-axis): Stack the two 4x4x2 pieces to form a 4x4x4 block. Cut at the 2cm mark.
    stage1_y_cuts = passes_per_4cm_cut
    
    # Cut 3 (x-axis): Stack the four 4x2x2 pieces to form a 4x4x4 block. Cut at the 2cm mark.
    stage1_x_cuts = passes_per_4cm_cut
    
    stage1_total_cuts = stage1_z_cuts + stage1_y_cuts + stage1_x_cuts
    
    print("### Stage 1: 4x4x4 -> 2x2x2 cubes ###")
    print(f"To make a central cut on a 4cm high block requires {passes_per_4cm_cut} passes.")
    print(f"Cuts for the central z-plane: {stage1_z_cuts}")
    print(f"Cuts for the central y-plane: {stage1_y_cuts}")
    print(f"Cuts for the central x-plane: {stage1_x_cuts}")
    print(f"Total cuts for Stage 1 = {stage1_z_cuts} + {stage1_y_cuts} + {stage1_x_cuts} = {stage1_total_cuts}\n")
    
    # --- Stage 2: Cut the eight 2x2x2 cubes into sixty-four 1x1x1 cubes ---
    # This requires cutting each 2cm dimension in half. For each axis, this corresponds to two planes of cuts (e.g., at x=1 and x=3).
    # We can stack the pieces to make a 4cm high block for each set of cuts.
    
    # X-axis cuts: Arrange the 8 cubes (2x2x2) into a 4x4x4 block. Make two cuts.
    stage2_x_cuts = passes_per_4cm_cut + passes_per_4cm_cut
    
    # Y-axis cuts: Rearrange the 16 pieces (1x2x2) into a 4x4x4 block. Make two cuts.
    stage2_y_cuts = passes_per_4cm_cut + passes_per_4cm_cut
    
    # Z-axis cuts: Rearrange the 32 pieces (1x1x2) into a 4x4x4 block. Make two cuts.
    stage2_z_cuts = passes_per_4cm_cut + passes_per_4cm_cut
    
    stage2_total_cuts = stage2_x_cuts + stage2_y_cuts + stage2_z_cuts
    
    print("### Stage 2: 2x2x2 -> 1x1x1 cubes ###")
    print("For each axis, we stack pieces into a 4cm high block and make two cuts.")
    print(f"Cuts for x-planes (at 1cm and 3cm marks): {passes_per_4cm_cut} + {passes_per_4cm_cut} = {stage2_x_cuts}")
    print(f"Cuts for y-planes (at 1cm and 3cm marks): {passes_per_4cm_cut} + {passes_per_4cm_cut} = {stage2_y_cuts}")
    print(f"Cuts for z-planes (at 1cm and 3cm marks): {passes_per_4cm_cut} + {passes_per_4cm_cut} = {stage2_z_cuts}")
    print(f"Total cuts for Stage 2 = {stage2_x_cuts} + {stage2_y_cuts} + {stage2_z_cuts} = {stage2_total_cuts}\n")

    # --- Total Calculation ---
    total_cuts = stage1_total_cuts + stage2_total_cuts
    
    print("### Total Minimum Cuts ###")
    print(f"Total = (Stage 1) + (Stage 2)")
    print(f"Total = ({stage1_z_cuts} + {stage1_y_cuts} + {stage1_x_cuts}) + ({stage2_x_cuts} + {stage2_y_cuts} + {stage2_z_cuts})")
    print(f"Total = {stage1_total_cuts} + {stage2_total_cuts} = {total_cuts}")

solve_cube_cutting()