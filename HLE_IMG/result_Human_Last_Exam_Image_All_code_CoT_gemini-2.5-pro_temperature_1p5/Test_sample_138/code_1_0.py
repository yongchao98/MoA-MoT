def solve_insect_sexing():
    """
    This function determines the indices for the biological sexes of the three insect pairs.
    - Pair A: Both insects are female Cuckoo bees with pointed abdomens. (F, F) -> Index 2
    - Pair B: The left wasp is a female (shorter antennae), and the right wasp is a male (yellow face, longer antennae). (F, M) -> Index 4
    - Pair C: The left bee is a male Long-horned bee (very long antennae), and the right is a female (shorter antennae). (M, F) -> Index 3
    """
    
    # Indices corresponding to the options: 1) M,M  2) F,F  3) M,F  4) F,M
    answer_A = 2  # F, F
    answer_B = 4  # F, M
    answer_C = 3  # M, F
    
    # Print the answer in the required format "A, B, C"
    print(f"{answer_A}, {answer_B}, {answer_C}")

solve_insect_sexing()