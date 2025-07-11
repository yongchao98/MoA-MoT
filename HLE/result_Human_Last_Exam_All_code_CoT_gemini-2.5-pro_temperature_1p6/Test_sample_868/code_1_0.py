def solve_dna_synthesis():
    """
    Calculates the composition of a newly synthesized DNA strand in a primer extension reaction.
    """
    # Step 1: Define the template and primer sequences.
    template = "GGACCGAATAACCGTAGAAGGCCGTAA"
    primer_3_to_5 = "TTGGCATCTTCC"

    print(f"Template strand: 5'-{template}-3'")
    print(f"Primer strand:   3'-{primer_3_to_5}-5'")
    print("-" * 50)

    # Step 2: Find the complement of the primer to locate the binding site.
    complement_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    primer_binding_seq = "".join([complement_map.get(base, 'N') for base in primer_3_to_5])
    
    print(f"The primer binds to its complementary sequence: 5'-{primer_binding_seq}-3'")
    
    # Locate the binding site on the template.
    binding_index = template.find(primer_binding_seq)

    if binding_index == -1:
        print("Error: Primer sequence does not bind to the template.")
        return

    print(f"The binding site is found at position {binding_index + 1} on the template.")
    print("-" * 50)
    
    # Step 3 & 4: Identify the template region that will be copied by the polymerase.
    # This is the region upstream (to the 5' end) of the primer binding site.
    template_to_copy = template[:binding_index]
    print(f"The polymerase extends from the 3' end of the primer, copying the upstream template region.")
    print(f"Template region to be copied: 5'-{template_to_copy}-3'")

    # Step 5: Determine the sequence of the newly synthesized DNA.
    synthesized_dna = "".join([complement_map.get(base, 'N') for base in template_to_copy])
    print(f"The newly synthesized (radiolabeled) DNA is: 3'-{synthesized_dna}-5'")
    print("-" * 50)

    # Step 6: Count the nucleotides in the newly synthesized DNA.
    count_A = synthesized_dna.count('A')
    count_T = synthesized_dna.count('T')
    count_C = synthesized_dna.count('C')
    count_G = synthesized_dna.count('G')

    # Final Result
    print("Final Answer:")
    print("The composition of the newly synthesized DNA is:")
    print(f"{count_A}A:{count_T}T:{count_C}C:{count_G}G")

solve_dna_synthesis()
<<<I>>>