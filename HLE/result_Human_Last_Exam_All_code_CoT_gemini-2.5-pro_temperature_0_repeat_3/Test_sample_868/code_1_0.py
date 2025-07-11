def solve_dna_synthesis():
    """
    Calculates the composition of a newly synthesized DNA strand in a primer extension reaction.
    """
    template_strand = "GGACCGAATAACCGTAGAAGGCCGTAA"
    primer_3_to_5 = "TTGGCATCTTCC"

    print("--- Step-by-Step Analysis ---")
    print(f"Template (5'->3'): {template_strand}")
    print(f"Primer   (3'->5'):   {' ' * 20}{primer_3_to_5}")

    # Step 1: Find the sequence on the template that the primer binds to.
    # This is the complement of the primer sequence.
    complement_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    
    # The sequence on the template is complementary to the primer.
    # e.g., 3'-T...C-5' on primer binds to 5'-A...G-3' on template.
    binding_site_on_template = "".join([complement_map[base] for base in primer_3_to_5])

    print(f"\nStep 1: The primer binds to its complementary sequence on the template.")
    print(f"The binding sequence is 5'-{binding_site_on_template}-3'.")

    # Step 2: Locate the binding site and identify the template for new synthesis.
    try:
        binding_start_index = template_strand.find(binding_site_on_template)
        if binding_start_index == -1:
            print("Error: Primer sequence not found on the template.")
            return

        # The polymerase synthesizes from the 3' end of the primer, copying the
        # template region that is 5' (upstream) to the binding site.
        template_for_synthesis = template_strand[:binding_start_index]
        
        print(f"\nStep 2: The polymerase will copy the template region upstream (5') of the binding site.")
        print(f"Template region to be copied: 5'-{template_for_synthesis}-3'")

    except ValueError:
        print("Error: Primer sequence not found on the template.")
        return

    # Step 3: Determine the sequence of the newly synthesized DNA strand.
    # The new strand is complementary to the template_for_synthesis.
    # The polymerase reads the template 3'->5' and synthesizes 5'->3'.
    # So, we read template_for_synthesis backwards and find the complement.
    newly_synthesized_strand = "".join([complement_map[base] for base in reversed(template_for_synthesis)])
    
    print(f"\nStep 3: The new DNA strand is synthesized.")
    print(f"Newly synthesized strand (5'->3'): {newly_synthesized_strand}")

    # Step 4: Count the nucleotides in the new strand.
    count_A = newly_synthesized_strand.count('A')
    count_T = newly_synthesized_strand.count('T')
    count_C = newly_synthesized_strand.count('C')
    count_G = newly_synthesized_strand.count('G')

    print("\n--- Final Composition ---")
    print("The composition of the newly synthesized DNA is:")
    print(f"{count_A}A:{count_T}T:{count_C}C:{count_G}G")

solve_dna_synthesis()
<<<I>>>