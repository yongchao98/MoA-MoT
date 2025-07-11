def solve_dna_synthesis():
    """
    Calculates the composition of a newly synthesized DNA strand based on a template and a primer.
    """
    # 1. Define the template and primer sequences.
    template = "GGACCGAATAACCGTAGAAGGCCGTAA"
    # The primer is given 3' -> 5', which is how it's written when aligned to the template.
    primer_3_to_5 = "TTGGCATCTTCC"

    print(f"Template (5'->3'): {template}")
    print(f"Primer   (3'->5'): {primer_3_to_5}")
    
    # 2. Find the primer binding site on the template.
    # We need the 5'->3' complement of the primer sequence.
    complements = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    
    # Generate the complementary sequence. Since the primer is 3'->5' and the template is 5'->3',
    # we just need to complement each base in order.
    binding_sequence = "".join([complements[base] for base in primer_3_to_5])
    
    print(f"\nSearching for primer's complementary sequence ({binding_sequence}) on the template...")

    # Find the starting position of the binding site.
    bind_index = template.find(binding_sequence)
    
    if bind_index == -1:
        print("Error: Primer does not bind to the template.")
        return

    print(f"Primer binds at index {bind_index}.")

    # 3. Identify the part of the template that will be copied.
    # This is the region upstream (5' end) of the primer binding site.
    template_to_copy = template[:bind_index]
    
    if not template_to_copy:
        print("No new DNA is synthesized as the primer binds at the very 5' end of the template.")
        return
        
    print(f"\nTemplate segment to be copied: 5'-{template_to_copy}-3'")

    # 4. Determine the sequence of the new DNA strand.
    # It will be the complement of the template segment.
    new_dna_strand = "".join([complements[base] for base in template_to_copy])
    
    print(f"Newly synthesized DNA strand (read 3'->5'): {new_dna_strand}")

    # 5. Count the nucleotides in the new strand.
    count_a = new_dna_strand.count('A')
    count_t = new_dna_strand.count('T')
    count_c = new_dna_strand.count('C')
    count_g = new_dna_strand.count('G')
    
    # 6. Print the final composition.
    print("\nComposition of the newly synthesized DNA:")
    # The final format requested: XA:YT:ZC:WG
    print(f"The composition of the newly synthesized DNA with radiolabeled nucleotides will be: {count_a}A:{count_t}T:{count_c}C:{count_g}G")

solve_dna_synthesis()
<<<I>>>