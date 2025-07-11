def solve_reaction_mechanism():
    """
    Identifies and describes the two pericyclic reactions in the given
    thermal transformation.
    """

    # --- Step 1: Electrocyclic Ring Opening ---
    pi_electrons_step1 = 4
    stereochemistry_step1 = "conrotatory"
    reaction_type_step1 = "electrocyclic ring opening"

    # --- Step 2: Electrocyclic Ring Closure ---
    pi_electrons_step2 = 6
    stereochemistry_step2 = "disrotatory"
    reaction_type_step2 = "electrocyclic ring closure"

    # --- Print the description of the overall transformation ---
    print("The thermal transformation involves two sequential pericyclic reactions:")
    print(f"1. A {pi_electrons_step1}π {stereochemistry_step1} {reaction_type_step1}.")
    print(f"2. A {pi_electrons_step2}π {stereochemistry_step2} {reaction_type_step2}.")

solve_reaction_mechanism()