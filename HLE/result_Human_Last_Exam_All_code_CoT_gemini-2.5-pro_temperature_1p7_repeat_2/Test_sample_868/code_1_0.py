import collections

def solve_dna_synthesis():
    """
    Solves the DNA synthesis problem by simulating primer binding and extension.
    """
    # Step 1: Define the input sequences from the problem.
    template_strand = "GGACCGAATAACCGTAGAAGGCCGTAA"
    primer_sequence_3_to_5 = "TTGGCATCTTCC"

    # Define the base-pairing rules for DNA.
    dna_complement_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

    # Step 2: Find the primer's binding site on the template.
    # The primer is 3'->5', so its complement on the template is 5'->3'.
    # We find the complement of the primer sequence to get the binding sequence.
    binding_site_sequence = "".join([dna_complement_map[base] for base in primer_sequence_3_to_5])

    # Find the starting index of the binding site on the template.
    try:
        binding_site_start_index = template_strand.index(binding_site_sequence)
    except ValueError:
        print("Error: Primer binding site not found on the template strand.")
        return

    # Step 3: Identify the part of the template that will be copied.
    # This is the sequence upstream (5'-ward) of the binding site.
    template_for_synthesis = template_strand[:binding_site_start_index]

    # Step 4: Determine the sequence of the newly synthesized DNA.
    # This will be the complement of the template segment being copied.
    # The synthesized strand will be antiparallel to the template segment.
    # Let's count from the 5' to 3' version of the new strand for consistency.
    newly_synthesized_dna_5_to_3 = "".join([dna_complement_map[base] for base in reversed(template_for_synthesis)])

    # Step 5: Count the composition of the newly synthesized DNA.
    composition = collections.Counter(newly_synthesized_dna_5_to_3)
    count_A = composition.get('A', 0)
    count_T = composition.get('T', 0)
    count_C = composition.get('C', 0)
    count_G = composition.get('G', 0)

    # Print the final result in the requested format.
    print("The composition of the newly synthesized DNA with radiolabeled nucleotides will be: " +
          str(count_A) + "A:" +
          str(count_T) + "T:" +
          str(count_C) + "C:" +
          str(count_G) + "G")

solve_dna_synthesis()
<<<I>>>