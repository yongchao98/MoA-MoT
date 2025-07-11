def solve():
    """
    Analyzes the band structures of three graphene nanoribbons and generates a concatenated classification string.
    """

    # --- Analysis of Ribbon 1 ---
    # Edge: Shows flat bands near the Fermi level, characteristic of Zigzag (Z).
    # Width: Counting bands at k=0 gives 7 bands above E=0. N=7.
    # Band: Bands are at the Fermi level (E=0), so it's metallic (0).
    ribbon1_edge = "Z"
    ribbon1_width = 7
    ribbon1_band = 0
    classification1 = f"{ribbon1_edge}{ribbon1_width}{ribbon1_band}"

    # --- Analysis of Ribbon 2 ---
    # Edge: No flat bands, gap closes at k=0, characteristic of Armchair (A).
    # Width: Counting bands at k!=0 gives 5 bands above E=0. N=5.
    # Band: Bands touch at the Fermi level, so it's metallic (0).
    ribbon2_edge = "A"
    ribbon2_width = 5
    ribbon2_band = 0
    classification2 = f"{ribbon2_edge}{ribbon2_width}{ribbon2_band}"

    # --- Analysis of Ribbon 3 ---
    # Edge: No flat bands, gap opens at k=0, characteristic of Armchair (A).
    # Width: Counting bands at k=0 gives 4 bands above E=0. N=4.
    # Band: Clear band gap exists, so it's semiconducting (1).
    ribbon3_edge = "A"
    ribbon3_width = 4
    ribbon3_band = 1
    classification3 = f"{ribbon3_edge}{ribbon3_width}{ribbon3_band}"

    # Concatenate the classifications as per the format "Edge_Width_Band"
    final_string = classification1 + classification2 + classification3
    
    print(final_string)

solve()