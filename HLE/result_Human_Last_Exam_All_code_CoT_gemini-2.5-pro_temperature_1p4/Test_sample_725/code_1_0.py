def solve_covalency():
    """
    Determines the relative covalency of CeF6(2-) and CeCl6(2-)
    based on provided orbital overlap information.
    """
    # Step 1: State the fundamental principle of covalent bonding.
    principle = "Covalency, the sharing of electrons in a chemical bond, is directly proportional to the extent of orbital overlap. Greater overlap leads to stronger covalency."
    print("Step 1: Recall the relationship between orbital overlap and covalency.")
    print(f"Principle: {principle}\n")

    # Step 2: State the information given in the problem for each compound.
    overlap_cef6 = "greater"
    overlap_cecl6 = "lesser" # Relative to CeF6(2-)
    print("Step 2: Analyze the given information.")
    print(f"Orbital overlap in CeF6(2-): {overlap_cef6}")
    print(f"Orbital overlap in CeCl6(2-): {overlap_cecl6}\n")

    # Step 3: Apply the principle to the given information to draw a conclusion.
    final_answer = ""
    if overlap_cef6 == "greater":
        final_answer = "stronger"
    else:
        # This case is not met by the problem statement but is included for logical completeness.
        final_answer = "weaker"

    print("Step 3: Formulate the conclusion.")
    print(f"Because the orbital overlap in CeF6(2-) is {overlap_cef6} than in CeCl6(2-), the CeF6(2-) complex will display {final_answer} covalency.")

solve_covalency()