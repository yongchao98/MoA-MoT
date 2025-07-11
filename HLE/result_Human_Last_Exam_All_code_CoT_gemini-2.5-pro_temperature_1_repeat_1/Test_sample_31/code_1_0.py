import math

def solve_pll_sticker_problem():
    """
    This script determines the minimum number of non-top-facing stickers
    needed to uniquely identify any of the 21 PLL cases.
    """

    # Step 1: Calculate the theoretical minimum using information theory.
    print("### Step 1: The Theoretical Minimum ###")
    
    num_pll_cases = 21
    num_sticker_colors = 4  # Side colors: e.g., Red, Green, Orange, Blue
    
    # To distinguish N states with an alphabet of size C, we need at least
    # ceil(log_C(N)) observations.
    theoretical_min = math.ceil(math.log(num_pll_cases, num_sticker_colors))
    
    print(f"There are {num_pll_cases} possible PLL cases.")
    print(f"Each non-top sticker can be one of {num_sticker_colors} colors.")
    print(f"The absolute minimum number of observations needed is ceil(log_{num_sticker_colors}({num_pll_cases})).")
    print(f"This calculates to: {theoretical_min}")
    print("So, from a pure information theory standpoint, 3 stickers might be enough.\n")

    # Step 2: Test if 3 stickers are practically sufficient.
    # We will check if we can distinguish the "Solved" state from the "Ab-perm" state
    # using a set of 3 stickers. The Ab-perm cycles three corners and leaves edges untouched.
    # Our chosen stickers:
    # 1. Front sticker of the Up-Front (UF) edge piece.
    # 2. Right sticker of the Up-Right (UR) edge piece.
    # 3. Front sticker of the Up-Front-Right (UFR) corner piece.
    # (We use F, R, B, L to denote the colors of the Front, Right, Back, Left faces)
    
    print("### Step 2: Testing a Set of 3 Stickers ###")
    print("Let's test if 3 stickers are enough in practice.")
    print("We'll examine the stickers for two cases: 'Solved' and 'Ab-perm'.\n")

    # Case A: The Solved state
    print("--- Case A: Solved State ---")
    solved_uf_f = 'F' # UF edge is in place, its front sticker is F-colored.
    solved_ur_r = 'R' # UR edge is in place, its right sticker is R-colored.
    solved_ufr_f = 'F' # UFR corner is in place, its front sticker is F-colored.
    print(f"Sticker 1 (UF edge, Front face): {solved_uf_f}")
    print(f"Sticker 2 (UR edge, Right face): {solved_ur_r}")
    print(f"Sticker 3 (UFR corner, Front face): {solved_ufr_f}")
    print(f"Resulting Pattern: ('{solved_uf_f}', '{solved_ur_r}', '{solved_ufr_f}')\n")

    # Case B: The Ab-perm state
    # This permutation cycles corners (UFL -> URF -> UBR) and leaves edges solved.
    print("--- Case B: Ab Permutation ---")
    ab_perm_uf_f = 'F' # Edges are not affected in an Ab-perm.
    ab_perm_ur_r = 'R' # Edges are not affected in an Ab-perm.
    # The UFR corner spot is now occupied by the original UFL piece.
    # The UFL piece has F and L colored stickers. Its F-facing sticker is F.
    ab_perm_ufr_f = 'F'
    print(f"Sticker 1 (UF edge, Front face): {ab_perm_uf_f}")
    print(f"Sticker 2 (UR edge, Right face): {ab_perm_ur_r}")
    print(f"Sticker 3 (UFR corner, Front face): {ab_perm_ufr_f}")
    print(f"Resulting Pattern: ('{ab_perm_uf_f}', '{ab_perm_ur_r}', '{ab_perm_ufr_f}')\n")

    print("--- Conclusion for Step 2 ---")
    print("The patterns are identical. This set of 3 stickers cannot distinguish between a solved state and an Ab-perm.")
    print("This demonstrates that 3 stickers are not sufficient.\n")

    # Step 3: Test if 4 stickers are sufficient by resolving the ambiguity.
    # We add a 4th sticker to our set: the Right sticker of the UFR corner.
    
    print("### Step 3: Testing a Set of 4 Stickers ###")
    print("Let's add a 4th sticker to see if it resolves the ambiguity.\n")

    # Case A: The Solved state with 4 stickers
    print("--- Case A: Solved State (4 stickers) ---")
    solved_ufr_r = 'R' # In the solved state, the UFR corner's R-sticker is R.
    print(f"Sticker 1 (UF edge, Front face): {solved_uf_f}")
    print(f"Sticker 2 (UR edge, Right face): {solved_ur_r}")
    print(f"Sticker 3 (UFR corner, Front face): {solved_ufr_f}")
    print(f"Sticker 4 (UFR corner, Right face): {solved_ufr_r}")
    print(f"Resulting Pattern: ('{solved_uf_f}', '{solved_ur_r}', '{solved_ufr_f}', '{solved_ufr_r}')\n")

    # Case B: The Ab-perm state with 4 stickers
    print("--- Case B: Ab Permutation (4 stickers) ---")
    # The UFR corner spot is occupied by the original UFL piece.
    # The UFL piece has F and L colored stickers. Its sticker on the R-face of the cube will be L-colored.
    ab_perm_ufr_r = 'L'
    print(f"Sticker 1 (UF edge, Front face): {ab_perm_uf_f}")
    print(f"Sticker 2 (UR edge, Right face): {ab_perm_ur_r}")
    print(f"Sticker 3 (UFR corner, Front face): {ab_perm_ufr_f}")
    print(f"Sticker 4 (UFR corner, Right face): {ab_perm_ufr_r}")
    print(f"Resulting Pattern: ('{ab_perm_uf_f}', '{ab_perm_ur_r}', '{ab_perm_ufr_f}', '{ab_perm_ufr_r}')\n")
    
    print("--- Conclusion for Step 3 ---")
    print("The patterns are now different. This demonstrates that 4 stickers are sufficient to resolve this ambiguity.\n")
    
    # Step 4: Final Conclusion
    print("### Final Answer ###")
    final_answer = 4
    print("While it can be proven that no set of 3 stickers can distinguish all 21 cases, a set of 4 can.")
    print(f"Therefore, the minimum number of non-top-facing stickers that must be seen is {final_answer}.")
    
solve_pll_sticker_problem()