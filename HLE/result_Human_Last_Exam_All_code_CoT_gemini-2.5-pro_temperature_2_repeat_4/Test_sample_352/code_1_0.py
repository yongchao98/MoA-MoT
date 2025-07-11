def design_mutagenesis():
    """
    Analyzes the negatively charged patch in protein x and suggests
    mutations to neutralize it for a site-directed mutagenesis experiment.
    """
    # Define the original amino acids and their positions
    original_patch = {
        47: ('Serine', 'S'),
        48: ('Glutamate', 'E'),
        49: ('Glutamate', 'E'),
        50: ('Aspartate', 'D')
    }

    # Define the ideal replacement amino acid
    replacement_aa = ('Alanine', 'A')

    # Explain the rationale
    print("--- Mutagenesis Design Plan ---")
    print("Goal: To relieve the inhibitory effect of the negatively charged patch at positions 47-50.\n")
    print("Rationale:")
    print("The original sequence S-E-E-D is highly negatively charged, especially with Serine at position 47 being a phosphorylation site.")
    print(f"To test the hypothesis that this charge causes autoinhibition, we will replace each amino acid with {replacement_aa[0]} ({replacement_aa[1]}).")
    print("Alanine is chosen because it is small, chemically inert, and uncharged, which will remove the negative charge with minimal structural disruption.\n")

    # Print the specific mutations
    print("--- Proposed Mutations ---")
    original_sequence = []
    mutant_sequence = []
    for position, (name, code) in original_patch.items():
        print(f"Position {position}: Replace {name} ({code}) -> {replacement_aa[0]} ({replacement_aa[1]})")
        original_sequence.append(code)
        mutant_sequence.append(replacement_aa[1])

    # Print the final sequences for clarity
    print("\n--- Final Sequence Change ---")
    original_str = "".join(original_sequence)
    mutant_str = "".join(mutant_sequence)
    print(f"Original Segment (47-50): {original_str}")
    print(f"Proposed Mutant Segment (47-50): {mutant_str}")

# Execute the function to get the plan and result
design_mutagenesis()
