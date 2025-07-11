def solve_dna_synthesis():
    """
    Solves the primer extension problem by finding the binding site,
    determining the synthesized sequence, and counting its nucleotide composition.
    """
    # 1. Define the template and primer sequences from the problem.
    template_5_3 = "GGACCGAATAACCGTAGAAGGCCGTAA"
    primer_3_5_seq = "TTGGCATCTTCC"

    # Helper function to find the complementary DNA sequence.
    def get_complement(dna_sequence):
        complement_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        return "".join([complement_map.get(base, '') for base in dna_sequence])

    # 2. Find the primer binding site.
    # The template's complement is generated in the 3'->5' direction.
    template_complement_3_5 = get_complement(template_5_3)
    
    # The primer is also 3'->5', so we can search for its sequence directly.
    try:
        binding_start_index = template_complement_3_5.index(primer_3_5_seq)
    except ValueError:
        print("Error: Primer sequence not found on the template's complement.")
        return

    # 3. Determine the part of the template that will be copied.
    # This is the section of the template from its 5' end up to where the primer binds.
    template_to_be_copied = template_5_3[0:binding_start_index]

    # 4. The newly synthesized DNA is complementary to this copied template region.
    newly_synthesized_dna_3_5 = get_complement(template_to_be_copied)
    
    # For standard representation and counting, let's use the 5'->3' sequence.
    newly_synthesized_dna_5_3 = newly_synthesized_dna_3_5[::-1]

    # 5. Count the nucleotides in the new DNA.
    count_A = newly_synthesized_dna_5_3.count('A')
    count_T = newly_synthesized_dna_5_3.count('T')
    count_C = newly_synthesized_dna_5_3.count('C')
    count_G = newly_synthesized_dna_5_3.count('G')

    # Output the results clearly.
    print(f"Template (5'->3'):         {template_5_3}")
    print(f"Primer (3'->5'):           {' ' * binding_start_index}{primer_3_5_seq}")
    print(f"Template part to be copied: 5'-{template_to_be_copied}-3'")
    print(f"Newly synthesized DNA:      5'-{newly_synthesized_dna_5_3}-3'")
    print("\n---")
    print("Final Answer:")
    print(f"The composition of the newly synthesized DNA with radiolabeled nucleotides will be: {count_A}A:{count_T}T:{count_C}C:{count_G}G")

solve_dna_synthesis()
<<<I>>>