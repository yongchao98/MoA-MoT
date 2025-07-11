def solve_primer_extension():
    """
    Calculates the composition of a newly synthesized DNA strand in a primer extension reaction.
    """
    # 1. Define the input sequences
    template_strand = "GGACCGAATAACCGTAGAAGGCCGTAA"
    primer_3_to_5 = "TTGGCATCTTCC"
    print(f"Template Strand (5'->3'): {template_strand}")
    print(f"Primer Strand (3'->5'):   {primer_3_to_5}")
    print("-" * 30)

    # 2. Determine the primer's complementary sequence to find its binding site
    complement_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    # To get the 5'->3' complementary sequence, we read the primer from right to left
    # and find the complement of each base.
    primer_complement_5_to_3 = "".join([complement_map[base] for base in reversed(primer_3_to_5)])

    # 3. Find the binding site on the template
    binding_site_index = template_strand.find(primer_complement_5_to_3)
    if binding_site_index == -1:
        print("Error: Primer does not bind to the template.")
        return

    print(f"The primer binds to the sequence: 5'-{primer_complement_5_to_3}-3'")

    # 4. Identify the portion of the template that will be copied by DNA polymerase
    # This is the region upstream (5'-ward) of the binding site.
    template_for_extension = template_strand[:binding_site_index]
    print(f"Template region for extension: 5'-{template_for_extension}-3'")

    # 5. Synthesize the new DNA strand (reverse complement of the extension template)
    complement_strand = "".join([complement_map[base] for base in template_for_extension])
    newly_synthesized_strand = complement_strand[::-1] # Reverse to get 5'->3' direction
    print(f"Newly synthesized strand: 5'-{newly_synthesized_strand}-3'")
    print("-" * 30)
    
    # 6. Count the nucleotides in the newly synthesized strand
    count_A = newly_synthesized_strand.count('A')
    count_T = newly_synthesized_strand.count('T')
    count_C = newly_synthesized_strand.count('C')
    count_G = newly_synthesized_strand.count('G')

    # 7. Print the final composition
    print("Final Composition of the Newly Synthesized DNA:")
    print(f"{count_A}A:{count_T}T:{count_C}C:{count_G}G")

solve_primer_extension()
<<<I>>>