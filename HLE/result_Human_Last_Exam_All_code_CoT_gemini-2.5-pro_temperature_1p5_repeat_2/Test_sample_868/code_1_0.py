def solve_dna_synthesis():
    """
    Solves a molecular biology problem involving primer extension to determine
    the composition of the newly synthesized DNA strand.
    """
    # Given sequences
    template = "GGACCGAATAACCGTAGAAGGCCGTAA"
    primer_3_to_5 = "TTGGCATCTTCC" # Note: This is given 3' -> 5'

    print("Step 1: Analyzing the components")
    print(f"Template (5'->3'): {template}")
    print(f"Primer (3'->5'):   {primer_3_to_5}\n")

    # Define base pairing rules
    complement_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

    # Step 2: Find the primer binding site on the template.
    # The primer (3'->5') binds to its complement (5'->3') on the template.
    # We create the complementary sequence that would be found on the template.
    binding_site_sequence = "".join([complement_map[base] for base in primer_3_to_5])
    
    print("Step 2: Finding the primer binding site")
    print(f"The 5'->3' sequence on the template that the primer binds to is the complement of the primer.")
    print(f"Calculated binding sequence: 5'-{binding_site_sequence}-3'\n")

    # Find the location of the binding site on the template
    try:
        binding_index_start = template.find(binding_site_sequence)
        if binding_index_start == -1:
            print("Error: Primer does not bind to the template.")
            return
        
        binding_index_end = binding_index_start + len(binding_site_sequence)
        
        # Display the alignment
        print("Alignment of primer on template:")
        print(f"Template: {template[:binding_index_start]}-{template[binding_index_start:binding_index_end]}-{template[binding_index_end:]}")
        print(f"          {' ' * len(template[:binding_index_start])}|{' ' * len(binding_site_sequence)}|")
        # To show the primer alignment, we reverse its sequence to match the 5'->3' template direction
        primer_5_to_3 = primer_3_to_5[::-1]
        print(f"Primer:   {' ' * binding_index_start}{primer_5_to_3}\n")


    except KeyError:
        print("Error: Invalid character in one of the sequences.")
        return

    # Step 3: Identify the part of the template to be copied.
    # DNA Polymerase extends from the 3' end of the primer, copying the template
    # region upstream (towards the 5' end of the template).
    template_to_copy = template[:binding_index_start]
    print("Step 3: Identifying the region for new DNA synthesis")
    print(f"Synthesis will copy the template region upstream of the binding site.")
    print(f"Template region to be copied: 5'-{template_to_copy}-3'\n")

    # Step 4: Determine the sequence of the newly synthesized DNA.
    # The new strand is complementary to the template region being copied.
    newly_synthesized_dna = "".join([complement_map[base] for base in template_to_copy])
    print("Step 4: Determining the new DNA sequence")
    print(f"The newly synthesized DNA strand is: 5'-{newly_synthesized_dna}-3'\n")

    # Step 5: Count the nucleotides in the new strand.
    count_A = newly_synthesized_dna.count('A')
    count_T = newly_synthesized_dna.count('T')
    count_C = newly_synthesized_dna.count('C')
    count_G = newly_synthesized_dna.count('G')
    
    print("Step 5: Calculating the final nucleotide composition")
    print(f"The composition of the newly synthesized DNA is:")
    print(f"{count_A}A:{count_T}T:{count_C}C:{count_G}G")


solve_dna_synthesis()
<<<I>>>