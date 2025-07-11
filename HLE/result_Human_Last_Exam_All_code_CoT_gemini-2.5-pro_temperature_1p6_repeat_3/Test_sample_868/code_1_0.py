def solve_dna_synthesis():
    """
    Calculates the composition of a newly synthesized DNA strand based on a template and a primer.
    """
    # 1. Define the template and primer sequences from the problem description.
    template_5_to_3 = "GGACCGAATAACCGTAGAAGGCCGTAA"
    primer_3_to_5 = "TTGGCATCTTCC"

    print(f"Template (5'->3'): {template_5_to_3}")
    print(f"Primer (3'->5'):   {primer_3_to_5}\n")

    # 2. Determine the primer's binding sequence (its 5'->3' complement).
    complement_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    primer_binding_seq_5_to_3 = "".join([complement_map[base] for base in primer_3_to_5])
    
    # The binding sequence is read in the same direction as the template (5'->3')
    # so we don't need to reverse it for the search.
    print(f"The primer binds to the complementary sequence: 5'-{primer_binding_seq_5_to_3}-3'")

    # 3. Find the binding site on the template.
    try:
        binding_start_index = template_5_to_3.index(primer_binding_seq_5_to_3)
        print(f"This sequence is found on the template starting at index {binding_start_index}.\n")
    except ValueError:
        print("Error: The primer does not bind to the template.")
        return

    # 4. Identify the portion of the template that will be copied by DNA polymerase.
    # This is the part of the template upstream (to the 5' side) of the primer binding site.
    template_to_be_copied = template_5_to_3[:binding_start_index]
    print(f"DNA polymerase will copy the template section: 5'-{template_to_be_copied}-3'")

    # 5. Synthesize the new strand by finding the reverse complement of the copied template section.
    # a. Reverse the sequence to simulate the 3'->5' reading direction of the polymerase.
    template_read_3_to_5 = template_to_be_copied[::-1]
    # b. Find the complement to generate the new 5'->3' strand.
    new_strand_5_to_3 = "".join([complement_map[base] for base in template_read_3_to_5])
    print(f"The newly synthesized DNA strand is: 5'-{new_strand_5_to_3}-3'\n")

    # 6. Count the nucleotides in the new strand.
    count_A = new_strand_5_to_3.count('A')
    count_T = new_strand_5_to_3.count('T')
    count_C = new_strand_5_to_3.count('C')
    count_G = new_strand_5_to_3.count('G')

    # 7. Print the final composition as an equation/ratio.
    print("The composition of the newly synthesized DNA is:")
    print(f"{count_A}A:{count_T}T:{count_C}C:{count_G}G")

solve_dna_synthesis()
<<<I>>>