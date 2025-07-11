def analyze_complex_lifetimes():
    """
    Analyzes the expected lifetimes of four Iridium complexes based on their structure.
    """
    # Step 1: Define the complexes and their key structural features.
    # The key feature for stability is the presence of an ortho-fluorine
    # on the phenyl ring of the cyclometalating ligand.
    complex_data = {
        1: {'name': 'Complex 1', 'has_ortho_fluorine': True},
        2: {'name': 'Complex 2', 'has_ortho_fluorine': False},
        3: {'name': 'Complex 3', 'has_ortho_fluorine': True},
        4: {'name': 'Complex 4', 'has_ortho_fluorine': True},
    }

    print("Principle: The lifetime of these Iridium emitters is strongly influenced by their chemical stability.")
    print("A fluorine atom at the ortho-position (C2) of the phenyl ring of the phenylpyridine ligand is known to decrease stability, leading to shorter lifetimes.\n")

    # Step 2: Identify complexes with the destabilizing feature.
    shorter_lifetime_complexes = []
    print("Analysis:")
    for number, data in complex_data.items():
        if data['has_ortho_fluorine']:
            status = "present -> shorter lifetime expected."
            shorter_lifetime_complexes.append(number)
        else:
            status = "absent -> longer lifetime expected."
        print(f"- Complex {number}: Ortho-fluorine is {status}")

    # Step 3: Print the final conclusion.
    shorter_lifetime_complexes.sort()
    
    # Format the final list of numbers as requested
    final_numbers = ", ".join(map(str, shorter_lifetime_complexes))

    print("\nConclusion: The complexes expected to show shorter lifetimes are those with an ortho-fluorine.")
    print(f"The numbers of these complexes are: {final_numbers}")

# Run the analysis
analyze_complex_lifetimes()