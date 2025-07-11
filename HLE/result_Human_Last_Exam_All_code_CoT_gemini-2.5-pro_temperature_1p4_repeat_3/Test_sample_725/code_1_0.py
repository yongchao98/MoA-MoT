def solve_covalency_problem():
    """
    Analyzes the relationship between orbital overlap and covalency
    for CeF6(2-) and CeCl6(2-).
    """
    # Covalency is directly proportional to the extent of orbital overlap.
    # Greater overlap implies more significant sharing of electron density.
    overlap_CeF6 = "greater"
    overlap_CeCl6 = "lesser"

    # The problem states that CeF6(2-) has greater orbital overlap than CeCl6(2-).
    if overlap_CeF6 == "greater":
        covalency_CeF6 = "stronger"
    else:
        covalency_CeF6 = "weaker"

    # Print the reasoning and the conclusion.
    print("Reasoning:")
    print("1. Covalency is determined by the extent of orbital overlap between the central metal and the ligands.")
    print("2. Greater orbital overlap leads to more effective electron sharing and therefore, stronger covalency.")
    print(f"3. It is given that the orbital overlap in CeF6(2-) is greater than in CeCl6(2-).")
    print("\nConclusion:")
    print(f"Therefore, CeF6(2-) displays {covalency_CeF6} covalency compared to CeCl6(2-).")

solve_covalency_problem()