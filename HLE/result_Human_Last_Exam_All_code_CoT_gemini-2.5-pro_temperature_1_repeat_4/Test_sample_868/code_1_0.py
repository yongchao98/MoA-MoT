def solve_primer_extension():
    """
    Solves the primer extension problem to find the composition of the new DNA strand.
    """
    # 1. Define the initial DNA sequences
    template_strand = "GGACCGAATAACCGTAGAAGGCCGTAA"
    primer_3_to_5 = "TTGGCATCTTCC"

    print(f"Template (5'->3'): {template_strand}")
    print(f"Primer (3'->5'):   {primer_3_to_5}")
    print("-" * 30)

    # 2. Determine the primer's binding site on the template
    complement_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    
    # The binding site on the template is the 5'->3' complement of the 3'->5' primer
    binding_site = "".join([complement_map[base] for base in primer_3_to_5])
    print(f"Step 1: The primer binds to its complementary sequence.")
    print(f"         The 5'->3' binding site on the template is: {binding_site}")

    # 3. Locate the binding site on the template
    try:
        start_index_of_binding = template_strand.index(binding_site)
        print(f"Step 2: The binding site was found on the template.")
    except ValueError:
        print("Error: The primer does not bind to the template strand.")
        return

    # Visualize the binding
    padding = " " * start_index_of_binding
    print(f"         Template: 5'-{template_strand}-3'")
    print(f"                   {padding}||||||||||||")
    print(f"         Primer:   3'-{primer_3_to_5}-5'")
    print("-" * 30)

    # 4. Identify the portion of the template to be copied
    # Synthesis proceeds from the 3' end of the primer, copying the template
    # upstream (towards the 5' end) of the binding site.
    template_to_be_copied = template_strand[:start_index_of_binding]
    print(f"Step 3: DNA Polymerase synthesizes a new strand complementary to the template")
    print(f"         region upstream of the primer.")
    print(f"         Template segment to be copied: 5'-{template_to_be_copied}-3'")

    # 5. Synthesize the new DNA strand
    newly_synthesized_dna = "".join([complement_map[base] for base in template_to_be_copied])
    print(f"         Newly synthesized DNA:         3'-{newly_synthesized_dna}-5'")
    
    # Read in standard 5'->3' direction by reversing the string
    newly_synthesized_dna_5_to_3 = newly_synthesized_dna[::-1]
    print(f"         (Reading 5'->3'):              5'-{newly_synthesized_dna_5_to_3}-3'")
    print("-" * 30)


    # 6. Count the nucleotides in the new strand
    # We count the bases in the new strand, which is `newly_synthesized_dna`
    count_A = newly_synthesized_dna.count('A')
    count_T = newly_synthesized_dna.count('T')
    count_C = newly_synthesized_dna.count('C')
    count_G = newly_synthesized_dna.count('G')

    print("Step 4: Count the bases in the newly synthesized DNA.")
    print(f"The composition of the newly synthesized DNA ({newly_synthesized_dna_5_to_3}) is:")
    print(f"{count_A}A:{count_T}T:{count_C}C:{count_G}G")

solve_primer_extension()
<<<I>>>