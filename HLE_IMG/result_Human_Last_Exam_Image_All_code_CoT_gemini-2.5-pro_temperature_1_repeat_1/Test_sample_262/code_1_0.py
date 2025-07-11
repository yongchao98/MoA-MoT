def analyze_complex_stability():
    """
    Analyzes the stability of four Iridium complexes based on their fluorination patterns
    to predict which will have shorter lifetimes in LECs.
    """
    # Define the properties of each complex based on the image
    # We assign a stability score based on fluorination.
    # Base assumption: Fluorination increases stability.
    # Refined assumption: Ortho-fluorination (at 2' position) is particularly stabilizing.
    # Scoring: +2 for ortho-F, +1 for para-F.
    complexes = {
        1: {'name': 'Complex 1', 'fluorination': "2',4'-difluoro", 'has_ortho_F': True, 'num_F': 2, 'score': 2 + 1},
        2: {'name': 'Complex 2', 'fluorination': "4'-fluoro", 'has_ortho_F': False, 'num_F': 1, 'score': 1},
        3: {'name': 'Complex 3', 'fluorination': "2'-fluoro", 'has_ortho_F': True, 'num_F': 1, 'score': 2},
        4: {'name': 'Complex 4', 'fluorination': "none", 'has_ortho_F': False, 'num_F': 0, 'score': 0},
    }

    print("Step 1: Assign stability scores based on fluorination.")
    print("Scoring rule: +2 for ortho-F, +1 for other F atoms.")
    for i in sorted(complexes.keys()):
        print(f"  - {complexes[i]['name']} (Fluorination: {complexes[i]['fluorination']}): Score = {complexes[i]['score']}")

    # Stability is proportional to the score. Lifetime is proportional to stability.
    # We are looking for complexes with shorter lifetimes, i.e., lower stability scores.
    # The major stabilizing feature is the ortho-fluorine. Complexes lacking it are less stable.
    shorter_lifetime_complexes = []
    print("\nStep 2: Identify complexes with shorter lifetimes.")
    print("Reasoning: Complexes lacking the highly stabilizing ortho-fluorine (2'-F) substituent will be less stable.")

    for i in sorted(complexes.keys()):
        if not complexes[i]['has_ortho_F']:
            shorter_lifetime_complexes.append(i)
            print(f"  - {complexes[i]['name']} does not have an ortho-F. It is expected to have a shorter lifetime.")
        else:
            print(f"  - {complexes[i]['name']} has an ortho-F. It is expected to have a longer lifetime.")

    print("\nConclusion:")
    print(f"The complexes expected to show shorter lifetimes are {shorter_lifetime_complexes}.")

analyze_complex_stability()