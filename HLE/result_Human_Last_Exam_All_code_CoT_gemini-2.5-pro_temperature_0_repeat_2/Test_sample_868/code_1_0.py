def solve_dna_synthesis():
    """
    This function calculates the composition of a newly synthesized DNA strand
    in a primer extension reaction.
    """
    # Step 1: Define the template and primer sequences
    template_strand = "GGACCGAATAACCGTAGAAGGCCGTAA"
    # The primer is given 3'->5', which is how it aligns for finding the binding site
    primer_3_to_5 = "TTGGCATCTTCC"

    print(f"Template (5'->3'): {template_strand}")
    print(f"Primer (3'->5'):   {primer_3_to_5}")
    print("-" * 30)

    # Step 2: Find the complementary sequence of the primer to locate the binding site
    complement_map = str.maketrans("ATCG", "TAGC")
    primer_binding_seq = primer_3_to_5.translate(complement_map)
    print(f"The primer binds to the complementary sequence (5'->3'): {primer_binding_seq}")

    # Step 3: Find the starting position of the binding site on the template
    try:
        binding_start_index = template_strand.find(primer_binding_seq)
        if binding_start_index == -1:
            print("Error: Primer does not bind to the template.")
            return
    except ValueError:
        print("Error: Primer does not bind to the template.")
        return

    # The part of the template that will be copied is upstream (5' side) of the binding site
    template_to_be_copied = template_strand[:binding_start_index]
    print(f"Template segment to be copied by polymerase: 5'-{template_to_be_copied}-3'")

    # Step 4: Synthesize the new strand by finding the complement of the copied template segment
    newly_synthesized_dna = template_to_be_copied.translate(complement_map)
    
    # To get the correct 5'->3' sequence, we need to reverse the result, as synthesis
    # proceeds by reading the template 3'->5'.
    # Template read direction (3'->5'): TAAGCCAGG
    # New strand synthesized (5'->3'): ATTCGGTCC
    newly_synthesized_dna_5_to_3 = newly_synthesized_dna[::-1]

    print(f"Newly synthesized DNA sequence (5'->3'): {newly_synthesized_dna_5_to_3}")
    print("-" * 30)

    # Step 5: Count the nucleotides in the newly synthesized strand
    count_A = newly_synthesized_dna_5_to_3.count('A')
    count_T = newly_synthesized_dna_5_to_3.count('T')
    count_C = newly_synthesized_dna_5_to_3.count('C')
    count_G = newly_synthesized_dna_5_to_3.count('G')

    # Step 6: Print the final composition
    print("Final Composition of Newly Synthesized DNA:")
    print(f"The composition of the newly synthesized DNA with radiolabeled nucleotides will be: {count_A}A:{count_T}T:{count_C}C:{count_G}G")

solve_dna_synthesis()
<<<I>>>