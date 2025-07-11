def analyze_complex_lifetimes():
    """
    Identifies which Iridium complexes are expected to have shorter lifetimes
    based on their chemical structure.
    """

    # Step 1: Define the structural information for each complex.
    # The key feature for stability is the position of fluorine atoms on the
    # cyclometalating phenylpyridine ligand, relative to the Ir-C bond.
    # We represent the positions as 'ortho', 'meta', or 'para'.
    complex_data = {
        1: {'fluorine_positions': ['meta', 'para']},
        2: {'fluorine_positions': ['ortho', 'para']},
        3: {'fluorine_positions': ['ortho']},
        4: {'fluorine_positions': ['ortho', 'meta', 'para']}
    }

    # Step 2: Explain the scientific principle.
    print("--- Analysis of Emitter Lifetime ---")
    print("The operational lifetime of these Iridium(III) emitters in LECs is largely")
    print("determined by their chemical stability.")
    print("\nA critical principle is that a fluorine (F) atom at the 'ortho' position")
    print("to the Iridium-Carbon (Ir-C) bond creates a known degradation pathway.")
    print("This leads to reduced stability and a shorter device lifetime.\n")
    print("We will now check each complex for the presence of an 'ortho'-fluorine.\n")

    # Step 3: Analyze each complex and identify those with shorter lifetimes.
    shorter_lifetime_complexes = []
    for complex_id, data in complex_data.items():
        has_ortho_fluorine = 'ortho' in data['fluorine_positions']
        if has_ortho_fluorine:
            shorter_lifetime_complexes.append(complex_id)
            status = "Shorter lifetime expected."
        else:
            status = "Longer lifetime expected."
        
        print(f"Complex {complex_id}: Has F at {data['fluorine_positions']} positions. -> {status}")

    # Step 4: Output the final result in the requested format.
    shorter_lifetime_complexes.sort()
    
    print("\n--- Conclusion ---")
    print("The complexes expected to show shorter lifetimes are those containing ortho-fluorinated ligands.")
    
    # Print the final list of numbers as a formatted string.
    final_equation = f"Shorter Lifetime Set = {{{', '.join(map(str, shorter_lifetime_complexes))}}}"
    print(final_equation)
    print(f"Number 1: {shorter_lifetime_complexes[0]}")
    print(f"Number 2: {shorter_lifetime_complexes[1]}")
    print(f"Number 3: {shorter_lifetime_complexes[2]}")


# Run the analysis
if __name__ == "__main__":
    analyze_complex_lifetimes()