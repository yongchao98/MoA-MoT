def solve_task():
    """
    This function provides the classification for the three pairs of images (A, B, C).
    The classification is based on identifying whether the species is in its native or invasive range.

    The mapping is as follows:
    1) both native
    2) both invasive
    3) right invaded
    4) left invaded

    - Pair A: Left image appears native, right image shows an invasive monoculture. -> 3
    - Pair B: Left image shows an invasive monoculture, right image appears native. -> 4
    - Pair C: Both images likely show a widely naturalized/invasive species outside its native range. -> 2
    """
    
    # Indices for A, B, and C
    index_A = 3
    index_B = 4
    index_C = 2
    
    # Print the result in the required format
    print(f"{index_A}, {index_B}, {index_C}")

solve_task()