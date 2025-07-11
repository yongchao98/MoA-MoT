def design_mutagenesis_experiment():
    """
    Prints the recommended amino acid replacements for the site-directed
    mutagenesis experiment on protein x.
    """
    # Define the target region and the proposed mutations
    mutations = {
        47: {'original': 'S (Serine)', 'mutant': 'A (Alanine)'},
        48: {'original': 'E (Glutamate)', 'mutant': 'A (Alanine)'},
        49: {'original': 'E (Glutamate)', 'mutant': 'A (Alanine)'},
        50: {'original': 'D (Aspartate)', 'mutant': 'A (Alanine)'}
    }

    print("To relieve the inhibitory effect of the negatively charged patch (aa 47-50), the recommended approach is to replace the native amino acids with Alanine.")
    print("\nThis strategy addresses the two sources of negative charge:")
    print("1. Mutation of Serine-47 to Alanine (S47A) prevents phosphorylation, thereby eliminating a potential negative charge at this site.")
    print("2. Mutation of Glutamate-48, Glutamate-49, and Aspartate-50 to Alanine (E48A, E49A, D50A) neutralizes their inherent negative charges.")

    print("\nThis creates a quadruple mutant (S47A/E48A/E49A/D50A), which is the most robust way to test the hypothesis that this specific patch's negative charge is autoinhibitory.")

    print("\nHere is the detailed plan for the mutagenesis:")
    
    # Using a f-string to create what could be considered the "final equation"
    # for the entire mutation from the original patch to the new one.
    original_sequence = "S(47)-E(48)-E(49)-D(50)"
    mutated_sequence = "A(47)-A(48)-A(49)-A(50)"

    print(f"Final Equation of Change: {original_sequence} -> {mutated_sequence}\n")
    
    # Print the specific replacements as per the request to output each number
    for position, changes in mutations.items():
      print(f"At position {position}, replace the original {changes['original']} with the mutant {changes['mutant']}.")

design_mutagenesis_experiment()