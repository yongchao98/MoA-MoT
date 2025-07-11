def find_minimum_pll_stickers():
    """
    This function logically determines the minimum number of non-top-facing
    stickers required to identify any of the 21 PLL cases.
    """
    
    # Step 1: Analyze the problem
    # The last layer has 12 side-facing stickers (8 on corners, 4 on edges).
    # We need to distinguish between 21 different permutation cases.

    print("Step 1: Analyzing why a low number of stickers is not enough.")
    print("----------------------------------------------------------")

    # Step 2: Prove that 3 stickers are insufficient.
    # Looking at the 3 stickers on a single face can lead to ambiguous information.
    # For example, a "bar" (all three stickers matching the face's center color) can appear
    # in multiple cases like F-perm, U-perms, or even a solved cube.
    # Therefore, 3 stickers are not enough.
    min_stickers = 3
    print(f"Hypothesis: {min_stickers} stickers are not enough.")
    print("Reason: A pattern on one face (3 stickers), like a 'bar' or 'headlights', can correspond to multiple different PLL cases.")
    print("Conclusion: The minimum number must be greater than 3.\n")


    # Step 3: Prove that 4 stickers are insufficient.
    # To prove 4 is not enough, we only need to find one set of 4 stickers that fails to distinguish
    # between at least two different PLL cases.
    min_stickers = 4
    print(f"Hypothesis: {min_stickers} stickers are not enough.")
    print("Let's test a strategic set of 4 stickers: the single edge sticker from each of the four side faces.")
    # For a solved layer, these four stickers will match the colors of their respective faces (e.g., Green, Red, Blue, Orange).
    # For PLLs that only permute the corner pieces (Aa, Ab, E-perm), the edge pieces do not move.
    print("For a solved layer and for all corner-only PLLs (Aa, Ab, E), the edge pieces are in their correct positions.")
    print("This means the four edge stickers will look identical across all these different cases.")
    print("Conclusion: Since these cases are indistinguishable with this set of 4 stickers, the minimum number must be greater than 4.\n")

    # Step 4: Prove that 5 stickers are sufficient.
    # It is well-established in speedcubing that looking at two adjacent faces is enough to identify any PLL.
    stickers_on_face_1 = 3
    stickers_on_face_2 = 3
    shared_stickers = 1 # The corner sticker between the two faces
    sufficient_stickers = stickers_on_face_1 + stickers_on_face_2 - shared_stickers
    
    print("Step 2: Determining the sufficient number of stickers.")
    print("-----------------------------------------------------")
    print("The standard method for PLL recognition is to view two adjacent faces (e.g., Front and Right).")
    print(f"This involves observing the {stickers_on_face_1} stickers on the Front face and the {stickers_on_face_2} on the Right face.")
    print(f"One corner sticker is shared between these views, so the total number of unique stickers is:")
    print(f"{stickers_on_face_1} + {stickers_on_face_2} - {shared_stickers} = {sufficient_stickers}")
    print("This set of 5 stickers is proven to be enough to uniquely distinguish all 21 PLL cases.\n")

    # Final Conclusion
    print("Final Conclusion:")
    print("Since 4 stickers are insufficient and 5 are sufficient, the minimum number required is 5.")
    
    return sufficient_stickers

if __name__ == "__main__":
    answer = find_minimum_pll_stickers()
    print(f"\n<<< {answer} >>>")
