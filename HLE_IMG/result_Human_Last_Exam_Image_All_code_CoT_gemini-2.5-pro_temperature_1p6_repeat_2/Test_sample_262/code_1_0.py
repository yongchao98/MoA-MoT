def analyze_complex_stability():
    """
    Analyzes the stability of four Iridium complexes to predict their lifetime in LECs.

    The stability and lifetime of these complexes are primarily affected by steric hindrance
    around the Iridium metal center. Substituents at the ortho-position of the
    cyclometalating phenyl ring cause steric strain, which weakens the crucial Ir-C bond,
    leading to faster degradation and shorter device lifetime.
    """

    # Representing the key structural feature for each complex.
    # 'ortho_F' is True if a Fluorine atom is at the sterically hindered ortho-position.
    complex_features = {
        1: {'name': 'Complex 1', 'ortho_F': False},
        2: {'name': 'Complex 2', 'ortho_F': False},
        3: {'name': 'Complex 3', 'ortho_F': True},
        4: {'name': 'Complex 4', 'ortho_F': True}
    }

    print("Analyzing the stability of the Iridium complexes...")
    print("-" * 40)
    print("Principle: The presence of a substituent (like Fluorine) at the ortho-position")
    print("of the phenyl ring relative to the Ir-C bond causes steric hindrance,")
    print("destabilizing the complex and shortening its lifetime.")
    print("-" * 40)

    shorter_lifetime_complexes = []
    for complex_id, features in complex_features.items():
        if features['ortho_F']:
            print(f"- {features['name']} ({complex_id}) has an ortho-Fluorine. Expected to have a shorter lifetime.")
            shorter_lifetime_complexes.append(complex_id)
        else:
            print(f"- {features['name']} ({complex_id}) lacks an ortho-Fluorine. Expected to be more stable.")

    shorter_lifetime_complexes.sort()
    
    print("\nConclusion:")
    print("The complexes expected to show shorter lifetimes are those with destabilizing ortho-substituents.")
    
    # The user instruction "output each number in the final equation!" is interpreted
    # as clearly stating the numbers of the identified complexes.
    result_string = " and ".join(map(str, shorter_lifetime_complexes))
    print(f"The identified complexes are: [{result_string}]")

analyze_complex_stability()