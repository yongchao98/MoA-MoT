def solve_plant_invasion_puzzle():
    """
    This function provides the solution to the plant invasion puzzle based on ecological analysis.

    - A (Garlic Mustard): The right image shows a dense monoculture, a hallmark of invasion, while the left shows a more diverse native community. This is 'right invaded' (3).
    - B (Lupine): The left image shows a massive, dense population characteristic of an invasion. The right shows a sparse population, more consistent with a native range. This is 'left invaded' (4).
    - C (Papaya): As a widely cultivated and naturalized species, both images show papaya in human-altered landscapes, likely outside its original native range. This is 'both invasive' (2).
    """
    
    # Assigning the index for each image pair
    # 1) both native
    # 2) both invasive
    # 3) right invaded
    # 4) left invaded
    
    index_A = 3
    index_B = 4
    index_C = 2

    # Print the result in the requested format "A, B, C"
    print(f"{index_A}, {index_B}, {index_C}")

solve_plant_invasion_puzzle()