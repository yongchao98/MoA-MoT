def solve_covalency():
    """
    Analyzes the relationship between orbital overlap and covalency to answer the question.
    """
    # The strength of a covalent bond is directly proportional to the extent of orbital overlap.
    # Greater overlap means more electron sharing and thus stronger covalency.
    overlap_CeF6_vs_CeCl6 = "greater"  # as given in the problem statement
    
    # Based on this relationship:
    if overlap_CeF6_vs_CeCl6 == "greater":
        covalency_conclusion = "stronger"
    else:
        covalency_conclusion = "weaker"

    print(f"Given that CeF6²⁻ has greater 4f orbital overlap with ligand orbitals than CeCl6²⁻, and because covalency is directly proportional to orbital overlap, CeF6²⁻ will display {covalency_conclusion} covalency.")

solve_covalency()