import itertools

def solve_pll_sticker_problem():
    """
    This function determines the minimum number of non-top-facing stickers
    required to uniquely identify any of the 21 PLL cases.
    """

    # --- Data Setup ---
    # The 12 non-top-facing stickers of the last layer.
    # We assign an index to each sticker position for reference.
    # cUFR_F = 0, cUFR_R = 1,  cUFL_F = 2, cUFL_L = 3,
    # cULB_L = 4, cULB_B = 5,  cUBR_B = 6, cUBR_R = 7,
    # eUF_F  = 8, eUL_L  = 9, eUB_B  = 10, eUR_R  = 11
    # F=0(Green), R=1(Red), L=2(Orange), B=3(Blue) - shifted for python
    # F=0(Green), R=1(Red), B=2(Blue), L=3(Orange)
    sticker_names = [
        "UFR_F", "UFR_R", "UFL_F", "UFL_L", "ULB_L", "ULB_B",
        "UBR_B", "UBR_R", "UF_F", "UL_L", "UB_B", "UR_R"
    ]

    # Pre-computed sticker colors for each of the 21 PLL cases.
    # Each row represents a PLL, and each column corresponds to a sticker's color.
    # This data is derived from the geometric transformations of each PLL permutation.
    PLL_DATA = {
        'Solved': [0, 1, 0, 3, 3, 2, 2, 1, 0, 3, 2, 1], # Solved
        'Aa':     [0, 1, 0, 3, 3, 2, 2, 1, 1, 0, 2, 3], # Edge Cycle (UF, UR, UL)
        'Ab':     [0, 1, 0, 3, 3, 2, 2, 1, 3, 1, 2, 0], # Edge Cycle (UF, UL, UR)
        'E':      [3, 0, 1, 2, 2, 1, 0, 3, 0, 3, 2, 1], # Corner Cycle (UFR-UFL-ULB-UBR)
        'H':      [2, 3, 2, 1, 0, 1, 0, 3, 2, 1, 0, 3], # Corner(opp) and Edge(opp) swap
        'Z':      [1, 2, 3, 0, 1, 2, 3, 0, 3, 2, 1, 0], # Corner(adj) and Edge(adj) swap
        'Ua':     [3, 0, 0, 1, 3, 2, 2, 3, 2, 0, 3, 1], # 3-Corner, 3-Edge cycle
        'Ub':     [0, 3, 2, 3, 3, 2, 1, 0, 3, 2, 0, 1], # 3-Corner, 3-Edge cycle
        'T':      [3, 0, 0, 3, 3, 2, 2, 1, 0, 1, 2, 3], # Corner(adj) swap, Edge(diag) swap
        'F':      [2, 3, 0, 1, 3, 2, 0, 3, 1, 3, 2, 0], # Corner(adj) swap, Edge(adj) swap
        'Ja':     [0, 3, 0, 1, 3, 2, 2, 1, 3, 0, 2, 1], # Corner(adj) swap, Edge(adj) swap
        'Jb':     [3, 0, 0, 3, 3, 2, 2, 1, 1, 3, 2, 0], # Corner(adj) swap, Edge(adj) swap
        'Ra':     [0, 1, 2, 1, 3, 0, 1, 2, 1, 3, 2, 0], # 3-Corner cycle, Edge(adj) swap
        'Rb':     [0, 1, 0, 3, 1, 2, 2, 3, 3, 1, 2, 0], # 3-Corner cycle, Edge(adj) swap
        'V':      [2, 3, 0, 3, 3, 0, 0, 1, 0, 1, 2, 3], # Corner(diag) swap, Edge(adj) swap
        'Na':     [2, 1, 0, 3, 2, 3, 0, 1, 2, 3, 0, 1], # Double adjacent corner swap
        'Nb':     [0, 3, 2, 1, 0, 1, 2, 3, 0, 1, 2, 3], # Double adjacent corner swap
        'Ga':     [0, 1, 2, 3, 1, 2, 2, 1, 2, 0, 3, 1], # G perm
        'Gb':     [1, 2, 0, 3, 3, 0, 2, 1, 3, 0, 2, 1], # G perm
        'Gc':     [3, 0, 0, 1, 3, 2, 1, 2, 2, 3, 0, 1], # G perm
        'Gd':     [0, 3, 0, 1, 3, 2, 2, 1, 1, 2, 3, 0], # G perm
    }
    
    pll_names = list(PLL_DATA.keys())
    num_stickers = len(sticker_names)
    
    # --- Computational Search ---
    for k in range(1, num_stickers + 1):
        print(f"Checking if {k} sticker(s) are enough to identify all 21 PLLs...")
        
        # Get all combinations of k sticker indices
        sticker_combinations = itertools.combinations(range(num_stickers), k)
        
        found_working_combination = False
        
        for combo in sticker_combinations:
            signatures = set()
            # Generate a signature for each PLL using the current sticker combo
            for pll in pll_names:
                signature = tuple(PLL_DATA[pll][i] for i in combo)
                signatures.add(signature)
            
            # If we have 21 unique signatures, this combo works
            if len(signatures) == 21:
                found_working_combination = True
                
                print("\nSuccess!")
                print(f"The combination of {k} stickers below is sufficient:")
                final_combo_names = [sticker_names[i] for i in combo]
                print(f"Sticker Combination: {final_combo_names}")
                print(f"The minimum number of stickers is therefore {k}.")

                # We print the equation showing how each piece contributes to the final number
                corner_stickers_in_combo = sum(1 for i in combo if i < 8)
                edge_stickers_in_combo = k - corner_stickers_in_combo

                print("\nThis final answer can be represented as:")
                print(f"Number of corner stickers seen: {corner_stickers_in_combo}")
                print(f"Number of edge stickers seen: {edge_stickers_in_combo}")
                print(f"Total stickers to identify PLL: {corner_stickers_in_combo} + {edge_stickers_in_combo} = {k}")
                
                return k

        if not found_working_combination:
            print(f"No combination of {k} stickers could distinguish all 21 cases.\n")

    return -1 # Should not be reached

if __name__ == "__main__":
    solve_pll_sticker_problem()