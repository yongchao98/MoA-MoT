def solve_insect_sexing():
    """
    This function determines the indices corresponding to the biological sexes
    of the three pairs of insects shown in the image.

    Analysis:
    - Pair A: The left insect has long antennae (male), the right has shorter antennae (female). This is M, F, which is option 3.
    - Pair B: The left insect is a female wasp (straight antennae, stinger). The right insect is a male wasp (curled antennae, yellow face). This is F, M, which is option 4.
    - Pair C: The left insect is a male long-horned bee (very long antennae). The right insect is a female (shorter antennae). This is M, F, which is option 3.
    """
    
    # Indices for the options provided:
    # 1) M, M
    # 2) F, F
    # 3) M, F
    # 4) F, M
    
    index_A = 3  # Male, Female
    index_B = 4  # Female, Male
    index_C = 3  # Male, Female
    
    # Print the result in the specified format "A, B, C"
    print(f"{index_A}, {index_B}, {index_C}")

solve_insect_sexing()