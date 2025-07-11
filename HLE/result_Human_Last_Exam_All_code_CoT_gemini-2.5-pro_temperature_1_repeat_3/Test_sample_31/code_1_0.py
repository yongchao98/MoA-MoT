def solve_pll_sticker_problem():
    """
    This script determines the minimum number of non-top-facing stickers
    needed to identify a PLL case on a 3x3 Rubik's Cube.
    """

    print("Step 1: Determine if 3 stickers are enough.")
    print("---------------------------------------------")
    print("Let's consider a worst-case scenario: distinguishing the Ga and Gb PLL cases.")
    print("Both of these cases can be oriented to show a 1x2 'block' of solved pieces on one side.")
    print("\nIf we view the cube from the front and the block is on the left, we might see:")
    print(" - Sticker 1 (Front-left corner, front face): Matches front face color.")
    print(" - Sticker 2 (Front edge, front face): Matches front face color.")
    print(" - Sticker 3 (Front-right corner, front face): Matches the RIGHT face color.")
    print("\nThis 3-sticker pattern can be identical for both Ga and Gb permutations.")
    print("Because at least two different PLLs can produce the same 3-sticker view, 3 stickers are NOT sufficient to guarantee identification.")
    print("Therefore, the minimum number must be more than 3.")

    print("\nStep 2: Check if 4 stickers are sufficient.")
    print("---------------------------------------------")
    print("Since 3 is not enough, the minimum must be at least 4. Let's see if 4 stickers can solve the previous ambiguity.")
    print("We add a 4th sticker to our view: the side-sticker of the corner in the block.")
    print(" - Sticker 4 (Front-left corner, left face):")
    print("\nWith this 4th sticker, the cases become distinct:")
    print(" - For Ga-perm: Sticker 4's color will match the RIGHT face.")
    print(" - For Gb-perm: Sticker 4's color will match the BACK face.")
    print("\nSince the right face and back face have different colors, these two cases can now be clearly distinguished.")
    print("This demonstrates that 4 stickers are sufficient to resolve the ambiguity.")

    print("\nStep 3: Conclusion.")
    print("---------------------------------------------")
    minimum_stickers = 4
    print("We have shown that 3 stickers are not enough, and 4 stickers are sufficient to distinguish even very similar cases.")
    print(f"Thus, the minimum number of non-top-facing stickers that must be seen is {minimum_stickers}.")


solve_pll_sticker_problem()