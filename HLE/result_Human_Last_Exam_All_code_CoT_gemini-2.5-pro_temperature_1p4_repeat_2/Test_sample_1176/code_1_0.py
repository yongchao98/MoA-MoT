def solve_geology_task():
    """
    Analyzes 10 statements about Cordilleran geology and classifies each as
    Consensus (C) or Debated (D).
    """

    # C = Consensus, D = Debated
    # Each entry in the list corresponds to a statement in order.
    classifications = [
        "C",  # (1) Morrison formation is a foredeep deposit.
        "D",  # (2) Metamorphic core complexes from a slab window.
        "D",  # (3) "Nevadaplano" plateau existence.
        "D",  # (4) "Arizonaplano" from magmatic thickening.
        "C",  # (5) Laramide structures are inboard of Sevier.
        "C",  # (6) Sierra Nevada/Idaho batholiths from Cordilleran arc.
        "C",  # (7) Ignimbrite flare-up propagated southwest.
        "D",  # (8) Colorado Plateau uplift history.
        "D",  # (9) Cascades arc elevation history.
        "D",  # (10) Laramide caused by Shatsky conjugate subduction.
    ]

    # Join the list of characters into a single string for the final answer.
    final_answer = "".join(classifications)

    # Print the resulting string.
    print(final_answer)

solve_geology_task()