def solve_dna_synthesis():
    """
    Calculates the composition of a newly synthesized DNA strand in a primer extension reaction.
    """
    template = "GGACCGAATAACCGTAGAAGGCCGTAA"
    primer_3_to_5 = "TTGGCATCTTCC"

    # Step 1: Find the sequence on the template where the primer binds.
    # We do this by finding the complement of the primer sequence.
    complement_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    
    # To get the binding site sequence in the 5'-3' direction, we reverse the 3'-5' primer
    # and then find the complement of each base.
    # However, a simpler way is to complement the 3'-5' primer to get its 5'-3' partner.
    binding_site_list = [complement_map[base] for base in reversed(primer_3_to_5)]
    binding_site = "".join(binding_site_list)

    # Step 2: Find the part of the template that will be copied.
    # This is the sequence immediately following the primer binding site.
    try:
        binding_index = template.find(binding_site)
        if binding_index == -1:
            print("Error: Primer does not bind to the template.")
            return
        
        extension_start_index = binding_index + len(binding_site)
        template_to_copy = template[extension_start_index:]
        
    except ValueError:
        print("Error: Primer binding site not found on the template.")
        return

    # Step 3: Determine the sequence of the newly synthesized (templated) DNA.
    # This new strand is complementary to the template_to_copy.
    new_strand_templated = "".join([complement_map[base] for base in template_to_copy])

    # Step 4: Account for the non-templated 'A' added by Taq polymerase.
    # The mention of 72C suggests a thermostable polymerase like Taq, which is known
    # to add a single Adenine to the 3' end of the product.
    # The new strand is synthesized 5' -> 3', so the non-templated 'A' is added at the end.
    # To correctly determine the 3' end, we need the 5'->3' sequence.
    # The complement of 5'-CCGTAA-3' (template_to_copy) is 3'-GGCATT-5'.
    # In standard 5'-to-3' notation, this is 5'-TTACGG-3'.
    # So we add the 'A' to this sequence.
    new_strand_5_to_3 = new_strand_templated[::-1]
    final_new_strand = new_strand_5_to_3 + "A"

    # Step 5: Count the nucleotides in the final synthesized strand.
    count_A = final_new_strand.count('A')
    count_T = final_new_strand.count('T')
    count_C = final_new_strand.count('C')
    count_G = final_new_strand.count('G')

    print("The part of the template to be copied is: 5'-{}'-3".format(template_to_copy))
    print("The newly synthesized DNA (including non-templated 'A') is: 5'-{}'-3".format(final_new_strand))
    print("\nThe composition of the newly synthesized DNA is:")
    print(f"{count_A}A:{count_T}T:{count_C}C:{count_G}G")

solve_dna_synthesis()
<<<F>>>