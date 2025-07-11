import collections

def solve_dna_synthesis():
    """
    This function calculates the composition of a newly synthesized DNA strand
    based on a template and a primer in a primer extension reaction.
    """
    # Step 1: Define the DNA template and primer sequences
    template_seq = "GGACCGAATAACCGTAGAAGGCCGTAA"
    primer_seq_3_to_5 = "TTGGCATCTTCC"

    # Define the rules for base complementarity
    complements = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

    # Step 2: Determine the sequence on the template where the primer binds.
    # The primer binds to its complement.
    # Primer (3'->5'): T T G G C A T C T T C C
    # Complementary template sequence (5'->3'): A A C C G T A G A A G G
    binding_site_on_template = "".join([complements[base] for base in reversed(primer_seq_3_to_5)])

    # Step 3: Find the binding site on the template and identify the region to be copied.
    # The polymerase copies the template region 5' (upstream) to the binding site.
    binding_index = template_seq.find(binding_site_on_template)

    if binding_index == -1:
        print("Primer does not bind to the template.")
        return

    template_to_be_copied = template_seq[:binding_index]

    # Step 4: Determine the sequence of the newly synthesized DNA strand.
    # The new strand is the complement of the template region being copied.
    newly_synthesized_strand = "".join([complements[base] for base in template_to_be_copied])

    # Step 5: Count the nucleotide composition of the new strand.
    composition = collections.Counter(newly_synthesized_strand)
    count_A = composition['A']
    count_T = composition['T']
    count_C = composition['C']
    count_G = composition['G']

    # Print the result in the specified format
    print("The newly synthesized DNA strand has the sequence (5'->3'):", "".join(reversed(newly_synthesized_strand)))
    print("The composition of the newly synthesized DNA with radiolabeled nucleotides will be:")
    # The final equation string showing each number
    print(f"{count_A}A:{count_T}T:{count_C}C:{count_G}G")

solve_dna_synthesis()