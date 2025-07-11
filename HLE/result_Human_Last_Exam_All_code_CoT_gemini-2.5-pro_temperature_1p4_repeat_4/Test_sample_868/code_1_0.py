import collections

# Step 1: Define the initial DNA sequences from the problem
template_strand = "GGACCGAATAACCGTAGAAGGCCGTAA"
primer_strand_3_to_5 = "TTGGCATCTTCC"

print(f"Template Strand (5'->3'): {template_strand}")
print(f"Primer Strand   (3'->5'):   {primer_strand_3_to_5}")
print("-" * 60)

# Step 2: Determine the primer's binding site on the template strand.
# The primer binds to its complementary sequence. We will create this sequence.
complement_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

# Find the complement of the 3'->5' primer. This gives the 5'->3' sequence on the template.
binding_site_on_template = "".join([complement_map[base] for base in primer_strand_3_to_5])
print(f"The primer binds to the following sequence on the template:")
print(f"Binding site (5'->3'): {binding_site_on_template}")
print("-" * 60)

# Step 3: Locate the binding site and identify the region to be copied.
try:
    binding_index = template_strand.index(binding_site_on_template)
    
    # The polymerase copies the part of the template before the binding site.
    template_to_be_copied = template_strand[:binding_index]
    
    print("Alignment:")
    print(f"Template: 5'-{template_to_be_copied}-|{binding_site_on_template}|-{template_strand[binding_index + len(binding_site_on_template):]}-3'")
    print(f"Primer:                         3'-{primer_strand_3_to_5}-5'")
    print("-" * 60)
    print(f"The region of the template to be copied is: 5'-{template_to_be_copied}-3'")
    print("-" * 60)

    # Step 4: Synthesize the new DNA strand.
    # The polymerase reads the template in the 3'->5' direction and synthesizes the new strand in the 5'->3' direction.
    # To simulate this, we reverse the template segment and find its complement.
    template_read_3_to_5 = template_to_be_copied[::-1]
    newly_synthesized_dna_5_to_3 = "".join([complement_map[base] for base in template_read_3_to_5])

    print(f"The newly synthesized DNA strand is (5'->3'): {newly_synthesized_dna_5_to_3}")
    print("-" * 60)
    
    # Step 5: Count the composition of the new strand.
    composition = collections.Counter(newly_synthesized_dna_5_to_3)
    count_A = composition.get('A', 0)
    count_T = composition.get('T', 0)
    count_C = composition.get('C', 0)
    count_G = composition.get('G', 0)
    
    print("Final Answer:")
    # The problem asks for the composition to be printed in a specific format.
    print(f"The composition of the newly synthesized DNA with radiolabeled nucleotides will be: {count_A}A:{count_T}T:{count_C}C:{count_G}G")

except ValueError:
    print("Error: The primer sequence does not have a complementary binding site on the template strand.")
<<<I>>>