def solve_insect_sexing():
    """
    This function provides the solution to the insect sexing puzzle based on entomological analysis.
    The analysis is as follows:
    - Pair A (Coelioxys bees): Left is Female (pointed abdomen), Right is Male (blunt, spined abdomen). This is option 4 (F, M).
    - Pair B (Polistes wasps): Left is Female (straight antennae), Right is Male (curled antennae). This is option 4 (F, M).
    - Pair C (Long-horned bees): Left is Male (very long antennae), Right is Female (short antennae). This is option 3 (M, F).
    """
    
    # Indices corresponding to the options: 1) M,M  2) F,F  3) M,F  4) F,M
    index_A = 4
    index_B = 4
    index_C = 3
    
    # Print the final answer in the required format
    print(f"{index_A}, {index_B}, {index_C}")

solve_insect_sexing()