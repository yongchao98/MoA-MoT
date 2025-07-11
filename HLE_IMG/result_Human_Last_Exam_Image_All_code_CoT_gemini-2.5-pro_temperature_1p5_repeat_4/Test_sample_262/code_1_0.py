def analyze_complex_stability():
    """
    Analyzes the stability of iridium complexes based on ligand fluorination
    and identifies those expected to have shorter lifetimes in LECs.
    """
    print("Step 1: Understand the factor determining lifetime.")
    print("The operational lifetime of these iridium emitters in LECs is primarily governed by their chemical stability. A key strategy to enhance stability is the fluorination of the cyclometalating (C^N) ligands.\n")

    print("Step 2: Analyze each complex based on fluorination.")
    analysis = {
        'Complex 1': 'Contains two difluorinated ligands. Expected to be very stable.',
        'Complex 2': 'Contains one non-fluorinated ligand. This is a known weak point for stability.',
        'Complex 3': 'Contains two monofluorinated ligands. This is the lowest degree of fluorination among the homoleptic complexes (1, 3, 4).',
        'Complex 4': 'Contains two trifluorinated ligands. Expected to be the most stable.'
    }
    for complex, note in analysis.items():
        print(f"- {complex}: {note}")
    print()

    print("Step 3: Conclude which complexes have shorter lifetimes.")
    print("Based on the analysis, a higher degree of fluorination correlates with higher stability and longer lifetime. The presence of a non-fluorinated ligand (Complex 2) is particularly detrimental.")
    print("Therefore, Complex 2 and Complex 3 are the least stable in the group and are expected to show shorter lifetimes.\n")

    shorter_lifetime_complexes = [2, 3]

    print("Final Answer: The complexes expected to show shorter lifetimes are:")
    for number in shorter_lifetime_complexes:
        print(number)

analyze_complex_stability()