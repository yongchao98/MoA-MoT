import sys

def solve_dna_synthesis():
    """
    Solves for the composition of a newly synthesized DNA strand in a primer extension reaction.
    """
    # Step 1: Define the oligo sequences involved in the reaction.
    # The longer oligo acts as the template.
    template = "GGACCGAATAACCGTAGAAGGCCGTAA"
    # The shorter oligo acts as the primer. It is given in 3' to 5' orientation.
    primer_3_to_5 = "TTGGCATCTTCC"

    # Define the rules for base pairing.
    complement_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

    # Step 2: Determine the sequence on the template that the primer binds to.
    # This is the complement of the primer sequence, read in the 5' to 3' direction.
    primer_binding_sequence = "".join([complement_map.get(base, 'N') for base in primer_3_to_5])

    # Step 3: Find the location of the binding site on the template strand.
    try:
        binding_start_index = template.index(primer_binding_sequence)
    except ValueError:
        print("Error: The primer sequence does not bind to the provided template.", file=sys.stderr)
        return

    # Step 4: Identify the portion of the template that will be copied.
    # DNA polymerase synthesizes new DNA complementary to the region upstream
    # (towards the 5' end) of the primer binding site.
    template_to_copy = template[:binding_start_index]

    # Step 5: Determine the sequence of the newly synthesized DNA.
    # The polymerase reads the template from 3' to 5'. So, we reverse the
    # template segment to be copied to simulate this reading direction.
    template_read_order = template_to_copy[::-1]
    
    # The new strand is synthesized in the 5' to 3' direction and is complementary
    # to the sequence being read.
    newly_synthesized_dna = "".join([complement_map.get(base, 'N') for base in template_read_order])

    # Step 6: Count the number of each nucleotide in the new DNA strand.
    count_A = newly_synthesized_dna.count('A')
    count_T = newly_synthesized_dna.count('T')
    count_C = newly_synthesized_dna.count('C')
    count_G = newly_synthesized_dna.count('G')

    # Step 7: Print the final composition in the required format.
    print(f"The template strand is: 5'-{template}-3'")
    print(f"The primer strand is: 3'-{primer_3_to_5}-5'")
    print(f"The sequence on the template where the primer binds is: 5'-{primer_binding_sequence}-3'")
    print(f"The template region to be copied is: 5'-{template_to_copy}-3'")
    print(f"The sequence of the newly synthesized DNA is: 5'-{newly_synthesized_dna}-3'")
    print("\n---")
    print("The composition of the newly synthesized DNA with radiolabeled nucleotides will be:")
    print(f"{count_A}A:{count_T}T:{count_C}C:{count_G}G")


solve_dna_synthesis()
<<<I>>>