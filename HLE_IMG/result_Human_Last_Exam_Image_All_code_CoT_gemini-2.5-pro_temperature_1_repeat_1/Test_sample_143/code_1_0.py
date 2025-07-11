def solve_plant_invasion():
    """
    This function provides the classification for the three image pairs based on visual evidence of invasive species behavior.

    - Pair A (Garlic Mustard): The right image shows a dense monoculture, a hallmark of an invasion. This is 'right invaded' (3).
    - Pair B (Lupine): The left image shows a dense, aggressive stand typical of an invasion. This is 'left invaded' (4).
    - Pair C (Papaya): Neither image shows signs of invasion; both appear to be in a native or cultivated non-invasive context. This is 'both native' (1).
    """
    
    # Indices corresponding to the analysis for pairs A, B, and C
    answer_A = 3
    answer_B = 4
    answer_C = 1
    
    # Print the result in the format "A, B, C"
    print(f"{answer_A}, {answer_B}, {answer_C}")

solve_plant_invasion()