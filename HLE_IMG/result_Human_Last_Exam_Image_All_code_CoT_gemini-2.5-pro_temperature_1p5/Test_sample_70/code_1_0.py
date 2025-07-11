def solve_reaction():
    """
    Identifies and describes the two pericyclic reactions in the given
    thermal transformation.
    """
    reaction1_electrons = 4
    reaction1_type = "electrocyclic ring opening"
    
    reaction2_electrons = 8
    reaction2_type = "electrocyclic ring closure"
    
    print("The thermal transformation proceeds through two sequential pericyclic reactions:")
    print(f"1. A {reaction1_electrons}-pi electron {reaction1_type}.")
    print(f"2. An {reaction2_electrons}-pi electron {reaction2_type}.")

solve_reaction()