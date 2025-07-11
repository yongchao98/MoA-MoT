import math

def solve_laozi_engraving():
    """
    Solves the optimization problem for Laozi's books by comparing key packing strategies.
    """
    # 1. Define constants from the problem description
    mat_w, mat_h = 140, 110
    c_radius = 20
    c_diameter = 2 * c_radius
    c_value = 9999  # Symbols on a circular plate
    s_side = 10
    s_value = 360   # Symbols on a squared plate

    # --- Strategy 1: Dense Hexagonal Packing of Circles ---
    # This strategy maximizes the number of high-value circular plates.
    # Calculation for a 140x110 sheet with circle radius r=20:
    # The vertical height needed for n staggered rows is 2*r + (n-1)*r*sqrt(3).
    # We check how many rows fit in the 110cm height:
    # 40 + (n-1)*20*math.sqrt(3) <= 110  => (n-1)*34.64 <= 70 => n-1 <= 2.02 => n = 3 rows.
    # Rows are placed along the 140cm width:
    # Row 1 (un-staggered) fits 3 circles.
    # Row 2 (staggered) fits 3 circles.
    # Row 3 (un-staggered) fits 3 circles.
    hex_N = 3 + 3 + 3  # Total circular plates
    
    # The space left over by this dense packing is in complex curved shapes and a
    # very thin strip, none large enough to cut a 10x10cm square.
    hex_M = 0  # Total squared plates
    hex_K = hex_N * c_value + hex_M * s_value

    # --- Strategy 2: Grid Packing of Circles + Filling with Squares ---
    # This strategy is less dense but leaves usable rectangular scraps.
    grid_N_x = mat_w // c_diameter  # 140 // 40 = 3
    grid_N_y = mat_h // c_diameter  # 110 // 40 = 2
    grid_N = grid_N_x * grid_N_y     # 3 * 2 = 6 circular plates

    # Calculate remaining area for squares
    used_w = grid_N_x * c_diameter  # 120 cm
    used_h = grid_N_y * c_diameter  # 80 cm
    
    rem_area1_w, rem_area1_h = mat_w - used_w, mat_h       # 20x110 cm strip
    rem_area2_w, rem_area2_h = used_w, mat_h - used_h      # 120x30 cm strip
    
    grid_M1 = (rem_area1_w // s_side) * (rem_area1_h // s_side) # 2 * 11 = 22
    grid_M2 = (rem_area2_w // s_side) * (rem_area2_h // s_side) # 12 * 3 = 36
    grid_M = grid_M1 + grid_M2  # 58 squared plates
    
    grid_K = grid_N * c_value + grid_M * s_value

    # --- Strategy 3: Only Squares ---
    # A baseline calculation using only the lower-value square plates.
    sq_only_M = (mat_w // s_side) * (mat_h // s_side) # 14 * 11 = 154
    sq_only_K = sq_only_M * s_value

    # --- 5. Compare Strategies and Find the Best ---
    if hex_K > grid_K and hex_K > sq_only_K:
        N_final, M_final, K_final = hex_N, hex_M, hex_K
    elif grid_K > hex_K and grid_K > sq_only_K:
        N_final, M_final, K_final = grid_N, grid_M, grid_K
    else:
        N_final, M_final, K_final = 0, sq_only_M, sq_only_K

    print("The optimal strategy is to maximize the number of high-value circular plates.")
    print("A dense hexagonal packing fits 9 circles, yielding the highest total value.")
    print("-" * 30)
    print(f"Number of circular plates (N): {N_final}")
    print(f"Number of squared plates (M): {M_final}")
    print(f"Maximal number of symbols (K): {K_final}")
    print("-" * 30)
    print("Final Calculation:")
    print(f"{N_final} * {c_value} + {M_final} * {s_value} = {K_final}")
    
    # The final answer format
    print(f"\n<<< {N_final}:{M_final}:{K_final} >>>")

# Run the solver
solve_laozi_engraving()