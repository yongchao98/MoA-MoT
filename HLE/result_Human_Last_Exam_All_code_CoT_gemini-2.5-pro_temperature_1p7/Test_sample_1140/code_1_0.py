def plan_peptide_ligation(total_length, num_fragments):
    """
    Calculates a fragmentation plan for Native Chemical Ligation.

    Args:
        total_length (int): The total length of the target peptide.
        num_fragments (int): The number of fragments to break the peptide into.
    """
    if num_fragments <= 0:
        print("Number of fragments must be positive.")
        return

    # Calculate the size of each fragment
    base_length = total_length // num_fragments
    remainder = total_length % num_fragments
    
    fragment_lengths = []
    for i in range(num_fragments):
        length = base_length
        # Distribute the remainder among the first few fragments
        if i < remainder:
            length += 1
        fragment_lengths.append(length)

    # Reverse to have larger fragments listed first for clarity
    fragment_lengths.sort(reverse=True)
    
    print(f"A strategy for synthesizing a {total_length}aa peptide using Native Chemical Ligation is to:")
    print(f"1. Synthesize {num_fragments} smaller peptide fragments via SPPS.")
    print(f"2. Ligate (join) them together to form the final product.")
    print("\nA possible fragmentation plan results in the following equation:")

    # Build and print the equation showing how fragments sum to the total
    equation_parts = [str(f) for f in fragment_lengths]
    equation_str = " + ".join(equation_parts) + f" = {total_length}"
    print(equation_str)


# Main execution for the specific problem
peptide_total_length = 100
fragments_to_create = 3
plan_peptide_ligation(peptide_total_length, fragments_to_create)
