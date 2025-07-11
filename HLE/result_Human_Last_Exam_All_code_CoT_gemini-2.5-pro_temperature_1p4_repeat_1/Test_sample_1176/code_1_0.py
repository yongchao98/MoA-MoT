def solve_geology_quiz():
    """
    Analyzes 10 statements about North American Cordilleran geology
    and determines if each is a consensus (C) or debated (D) view.
    """

    # Each tuple contains the statement number and its classification ('C' or 'D').
    analysis = [
        (1, 'C'),  # Consensus: Morrison Formation as a foredeep deposit.
        (2, 'D'),  # Debated: Cause of metamorphic core complexes (slab window is one of several hypotheses).
        (3, 'C'),  # Consensus: Existence of a Cretaceous "Nevadaplano" is the dominant model.
        (4, 'D'),  # Debated: The "Arizonaplano" and its magmatic thickening origin are highly debated.
        (5, 'C'),  # Consensus: Laramide structures are inboard of Sevier structures.
        (6, 'C'),  # Consensus: Sierra Nevada/Idaho batholiths are remnants of the ancestral Cordilleran arc.
        (7, 'C'),  # Consensus: Southwestward propagation of Eocene volcanism (Yellowstone hotspot track).
        (8, 'D'),  # Debated: Uplift history of the Colorado Plateau is a classic, major debate.
        (9, 'D'),  # Debated: The uplift timing of the modern Cascades arc is complex and debated.
        (10, 'D') # Debated: Cause of Laramide deformation (Shatsky conjugate is a leading but not sole hypothesis).
    ]

    # Combine the classifications into a single continuous string.
    result_string = "".join([item[1] for item in analysis])

    # Print the final result string.
    print(result_string)

solve_geology_quiz()