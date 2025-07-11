import math

def solve_cutting_problem():
    """
    This script calculates the minimum number of cuts to divide a 4x4x4 cube
    into 1x1x1 cubes with a knife that has a limited cutting depth.
    """
    # --- Problem Parameters ---
    cube_dim = 4
    knife_depth = 2

    # --- Introduction ---
    print("To solve this, we calculate the minimum cuts for each of the 3 dimensions and sum them up.")
    print("We can use stacking to cut multiple pieces at once.\n")

    # --- Analysis for a single dimension ---
    print(f"Analyzing the cuts for a single {cube_dim}cm dimension with a {knife_depth}cm deep knife:")

    # To get 4 pieces from 1, we ideally need log2(4) = 2 stages of cutting.
    
    # Stage 1: Cut the initial 4cm piece(s) in half.
    height_stage1 = 4
    passes_stage1 = math.ceil(height_stage1 / knife_depth)
    print(f"- Stage 1: Halve the {height_stage1}cm pieces.")
    print(f"  The material is {height_stage1}cm thick, requiring {int(passes_stage1)} passes.")

    # Stage 2: Cut the resulting 2cm pieces in half.
    height_stage2 = 2
    passes_stage2 = math.ceil(height_stage2 / knife_depth)
    print(f"- Stage 2: Halve the resulting {height_stage2}cm pieces.")
    print(f"  By stacking, the material to cut is {height_stage2}cm thick, requiring {int(passes_stage2)} pass.")
    
    # Total for one dimension
    cuts_per_dimension = int(passes_stage1 + passes_stage2)
    print(f"\nTotal passes for one dimension = {int(passes_stage1)} + {int(passes_stage2)} = {cuts_per_dimension}.\n")

    # --- Final Calculation ---
    cuts_x = cuts_per_dimension
    cuts_y = cuts_per_dimension
    cuts_z = cuts_per_dimension
    total_cuts = cuts_x + cuts_y + cuts_z

    print("--- Final Answer ---")
    print("This process is repeated for all three dimensions (X, Y, Z).")
    print("The total minimum number of cuts is the sum:")
    print(f"{cuts_x} + {cuts_y} + {cuts_z} = {total_cuts}")

solve_cutting_problem()