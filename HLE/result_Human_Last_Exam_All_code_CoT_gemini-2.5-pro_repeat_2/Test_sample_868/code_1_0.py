def solve_dna_synthesis():
    """
    Solves the primer extension problem by finding the primer binding site,
    determining the newly synthesized sequence, and calculating its base composition.
    """
    # Step 1: Define the DNA sequences
    template_strand = "GGACCGAATAACCGTAGAAGGCCGTAA"
    # The primer is 3’-TTGGCATCTTCC-5’, its sequence is TTGGCATCTTCC
    primer_sequence_3_to_5 = "TTGGCATCTTCC"

    # Step 2: Find the complementary sequence on the template that the primer binds to
    complement_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    
    # To find where the 3'->5' primer binds on the 5'->3' template, we find the
    # complement of the primer sequence.
    primer_binding_site_seq = "".join([complement_map.get(base, 'N') for base in primer_sequence_3_to_5])

    # Step 3: Locate the binding site on the template
    binding_start_index = template_strand.find(primer_binding_site_seq)

    if binding_start_index == -1:
        print("Error: Primer does not bind to the template.")
        return

    # Step 4: Identify the portion of the template that will be copied by DNA polymerase.
    # This is the region 5' (upstream) of the primer binding site.
    template_to_copy = template_strand[:binding_start_index]

    # Step 5: Determine the sequence of the newly synthesized DNA strand.
    # It is the complement of the template_to_copy.
    newly_synthesized_dna = "".join([complement_map.get(base, 'N') for base in template_to_copy])

    # Step 6: Count the bases in the newly synthesized DNA.
    count_A = newly_synthesized_dna.count('A')
    count_T = newly_synthesized_dna.count('T')
    count_C = newly_synthesized_dna.count('C')
    count_G = newly_synthesized_dna.count('G')

    # Print the final result in the required format
    print(f"The newly synthesized DNA sequence is: 5'-{newly_synthesized_dna}-3'")
    print("The composition of the newly synthesized DNA with radiolabeled nucleotides will be: {A}A:{T}T:{C}C:{G}G".format(
        A=count_A, T=count_T, C=count_C, G=count_G
    ))

solve_dna_synthesis()
<<<I>>>