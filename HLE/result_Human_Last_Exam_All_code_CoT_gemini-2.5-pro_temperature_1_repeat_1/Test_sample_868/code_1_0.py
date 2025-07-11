def solve_dna_synthesis():
    """
    Solves the primer extension problem by simulating DNA synthesis.
    """
    # 1. Define the template and primer sequences
    template_strand = "GGACCGAATAACCGTAGAAGGCCGTAA"
    primer_strand_3_to_5 = "TTGGCATCTTCC"

    print(f"Template (5'->3'): {template_strand}")
    print(f"Primer (3'->5'):   {primer_strand_3_to_5}\n")

    # 2. Determine the sequence on the template where the primer binds
    complement_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    
    # The complement to the 3'->5' primer is the 5'->3' binding site
    binding_site_on_template = "".join([complement_map[base] for base in primer_strand_3_to_5])
    print(f"Step 1: Find the primer's binding site on the template.")
    print(f"The 5'->3' sequence complementary to the primer is: {binding_site_on_template}\n")

    # 3. Find the starting position of the binding site on the template
    try:
        binding_index = template_strand.index(binding_site_on_template)
        print(f"Step 2: Locate the binding site on the template strand.")
        print(f"Template: {template_strand}")
        # Create a visual alignment
        alignment_str = ' ' * binding_index + '|' * len(binding_site_on_template)
        print(f"          {alignment_str}")
        print(f"Binding:  {' ' * binding_index}{binding_site_on_template}\n")
    except ValueError:
        print("Error: Primer binding site not found on the template.")
        return

    # 4. Identify the portion of the template that gets copied by the polymerase
    # This is the region 5' (upstream) to the primer binding site
    template_to_be_copied = template_strand[:binding_index]
    print(f"Step 3: Identify the template region to be synthesized.")
    print(f"Synthesis starts from the primer's 3' end and copies the region upstream (5'-ward) of the binding site.")
    print(f"Template region to be copied: 5'-{template_to_be_copied}-3'\n")

    # 5. Determine the sequence of the newly synthesized DNA
    newly_synthesized_dna = "".join([complement_map[base] for base in template_to_be_copied])
    print(f"Step 4: Determine the sequence of the new DNA.")
    print(f"The newly synthesized DNA strand is: 3'-{newly_synthesized_dna}-5'\n")

    # 6. Count the nucleotides in the new strand
    count_A = newly_synthesized_dna.count('A')
    count_T = newly_synthesized_dna.count('T')
    count_C = newly_synthesized_dna.count('C')
    count_G = newly_synthesized_dna.count('G')

    print(f"Step 5: Count the nucleotides in the newly synthesized DNA ({newly_synthesized_dna}).")
    print(f"Number of A's: {count_A}")
    print(f"Number of T's: {count_T}")
    print(f"Number of C's: {count_C}")
    print(f"Number of G's: {count_G}\n")
    
    # Final Answer
    print("Final Answer:")
    print(f"The composition of the newly synthesized DNA with radiolabeled nucleotides will be: {count_A}A:{count_T}T:{count_C}C:{count_G}G")

# Run the simulation
solve_dna_synthesis()
<<<I>>>