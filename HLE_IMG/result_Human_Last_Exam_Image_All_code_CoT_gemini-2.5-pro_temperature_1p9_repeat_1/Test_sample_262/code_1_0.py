def find_shorter_lifetime_complexes():
    """
    Identifies which of the four iridium complexes are expected to have shorter
    lifetimes based on their chemical structure.
    """

    # Step 1 & 2: Define the scientific principle.
    # The operational lifetime of these complexes is primarily dictated by their chemical stability.
    # A major degradation pathway involves the breaking of the Iridium-Carbon (Ir-C) bond.
    print("Scientific Principle:")
    print("The stability of cyclometalated iridium complexes, and thus their operational lifetime in devices, is highly sensitive to substituents on the phenyl ring.")
    print("Specifically, a fluorine atom at the 'ortho' position (adjacent to the Ir-C bond) introduces steric and electronic strain that weakens the Ir-C bond.")
    print("A weaker Ir-C bond leads to faster molecular degradation and, consequently, a shorter device lifetime.")
    print("-" * 40)

    # Step 3: Evaluate each complex based on this principle.
    # We will identify which complexes have this destabilizing feature.
    complexes_data = {
        1: {'feature': 'ortho and para fluorine', 'has_ortho_F': True},
        2: {'feature': 'para fluorine', 'has_ortho_F': False},
        3: {'feature': 'ortho fluorine', 'has_ortho_F': True},
        4: {'feature': 'no fluorine', 'has_ortho_F': False}
    }

    shorter_lifetime_ids = []
    print("Evaluation of Complexes:")
    for cid, data in complexes_data.items():
        if data['has_ortho_F']:
            print(f"Complex {cid}: Contains an ortho-fluorine. It is expected to have a SHORTER lifetime.")
            shorter_lifetime_ids.append(cid)
        else:
            print(f"Complex {cid}: Lacks an ortho-fluorine. It is expected to have a LONGER lifetime.")
    print("-" * 40)

    # Step 4: Final Conclusion.
    print("Conclusion:")
    print("The complexes with ortho-fluorine substituents are expected to show shorter lifetimes.")
    print("These are complexes:")
    for complex_id in shorter_lifetime_ids:
        print(complex_id)

find_shorter_lifetime_complexes()
<<<F>>>