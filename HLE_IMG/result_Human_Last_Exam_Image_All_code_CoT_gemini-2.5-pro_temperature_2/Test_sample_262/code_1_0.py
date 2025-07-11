def solve_chemistry_stability():
    """
    Analyzes the stability of four Iridium complexes to determine which have shorter lifetimes.

    The stability of these complexes in LECs is influenced by two main factors related to the
    fluorine (F) substituents on the phenylpyridine ligands:
    1.  Electronic Stabilization: Fluorine atoms are electron-withdrawing, which strengthens the
        critical Ir-C bond, increasing stability. More F atoms (at non-hindered positions) is better.
    2.  Steric Hindrance: An F atom at the 'ortho' position (next to the Ir-C bond) causes
        steric strain, which weakens the Ir-C bond and greatly reduces stability.

    We can model this with a simple scoring system:
    - Each F atom at a non-ortho position adds +1 to the stability score.
    - Each F atom at an ortho position adds -5 to the stability score (a large penalty).
    - A complex with no F atoms has a score of 0.

    A lower score indicates lower stability and a shorter expected lifetime.
    """

    # Data describing each complex based on the provided image.
    # Each complex has two identical C^N ligands.
    complexes = [
        {'id': 1, 'name': 'Complex 1', 'non_ortho_F_per_ligand': 2, 'ortho_F_per_ligand': 0},
        {'id': 2, 'name': 'Complex 2', 'non_ortho_F_per_ligand': 1, 'ortho_F_per_ligand': 0},
        {'id': 3, 'name': 'Complex 3', 'non_ortho_F_per_ligand': 0, 'ortho_F_per_ligand': 1},
        {'id': 4, 'name': 'Complex 4', 'non_ortho_F_per_ligand': 0, 'ortho_F_per_ligand': 0}
    ]

    # Points for scoring
    STABILIZING_F_SCORE = 1
    DESTABILIZING_ORTHO_F_SCORE = -5

    print("Calculating stability scores for each complex:")

    results = []
    for c in complexes:
        # Each complex has 2 cyclometalating ligands
        num_ligands = 2
        non_ortho_F_total = c['non_ortho_F_per_ligand'] * num_ligands
        ortho_F_total = c['ortho_F_per_ligand'] * num_ligands

        # Calculate the score based on our model
        score = (non_ortho_F_total * STABILIZING_F_SCORE) + (ortho_F_total * DESTABILIZING_ORTHO_F_SCORE)
        c['score'] = score
        results.append(c)
        
        # Print the equation and result for each complex
        print(f"{c['name']}: ( {non_ortho_F_total} non-ortho F * {STABILIZING_F_SCORE} ) + ( {ortho_F_total} ortho F * {DESTABILIZING_ORTHO_F_SCORE} ) = {score}")


    # Sort complexes by stability score to find the least stable ones
    results.sort(key=lambda x: x['score'])

    # The question asks for those with shorter lifetimes, which corresponds to the lowest stability.
    # Based on the scores, Complex 3 (-10) and Complex 4 (0) are the least stable.
    shortest_lifetime_complexes = [c for c in results if c['score'] <= 0]
    
    ids = sorted([c['id'] for c in shortest_lifetime_complexes])

    print("\nConclusion:")
    print(f"The complexes are ranked by stability score: {[f'Complex {c["id"]} ({c["score"]})' for c in results]}.")
    print(f"The complexes with the lowest stability scores, and thus expected shorter lifetimes, are Complexes {ids[0]} and {ids[1]}.")

solve_chemistry_stability()
<<<J>>>