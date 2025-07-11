def solve_insect_sex():
    """
    Solves the insect sex identification puzzle.
    
    A: Left bee has a pointed abdomen (female Coelioxys ovipositor). Right has a spined abdomen tip (male). -> F, M (4)
    B: Left wasp has curled antennae tips (male). Right has straight antennae (female). -> M, F (3)
    C: Left bee has very long antennae (male long-horned bee). Right has short antennae (female). -> M, F (3)
    """
    
    # Indices corresponding to the analysis for pairs A, B, and C
    index_A = 4
    index_B = 3
    index_C = 3
    
    # Format and print the final answer as a comma-separated string
    print(f"{index_A}, {index_B}, {index_C}")

solve_insect_sex()