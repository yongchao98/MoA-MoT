def check_mutation_answer():
    """
    Checks the correctness of the LLM's answer by translating the gene sequences
    and identifying the mutation most likely to cause a loss of function.
    """
    genetic_code = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_STOP_', 'TAG':'_STOP_',
        'TGC':'C', 'TGT':'C', 'TGA':'_STOP_', 'TGG':'W',
    }

    sequences = {
        "Intact gene": "ATGTTTCTCGCTGGTACTTCTGTGGATGAACATATTTATTGTCGT",
        "Mutant 1":    "ATGTTCTACGCTGGTACTTCTGTGGATGAACATATTTATTGTCGC",
        "Mutant 2":    "ATGTTCTAAGCTGGTACTTCTGTGGATGAACATATTTATTGTCGC",
        "Mutant 3":    "ATGTTTTACGCTGGTGTCACTTCTGTGGATGAACATATTTATTGTCGT",
        "Mutant 4":    "ATGTTTTACGCTACTTCTGTGGATGAACATATTTATTGTCGT"
    }

    llm_answer = "D" # Corresponds to Mutant 2

    def translate(dna):
        """Translates a DNA sequence into a protein sequence."""
        protein = []
        for i in range(0, len(dna) - (len(dna) % 3), 3):
            codon = dna[i:i+3]
            amino_acid = genetic_code.get(codon, 'X') # 'X' for unknown
            if amino_acid == '_STOP_':
                protein.append(amino_acid)
                break
            protein.append(amino_acid)
        return protein

    results = {name: translate(seq) for name, seq in sequences.items()}

    # The goal is to find the mutation most likely to eliminate the compound,
    # which means creating a non-functional protein. A nonsense mutation (premature stop)
    # is the most effective way to achieve this.

    # Check Mutant 2
    mutant2_protein = results["Mutant 2"]
    if "_STOP_" not in mutant2_protein or len(mutant2_protein) > 3:
        return "Incorrect: The analysis of Mutant 2 is wrong. It should introduce a premature stop codon early in the sequence. The translated protein is: " + "-".join(mutant2_protein)

    # Check if any other mutant also has a premature stop codon
    for name, protein in results.items():
        if name in ["Intact gene", "Mutant 2"]:
            continue
        if "_STOP_" in protein:
            return f"Incorrect: The analysis is incomplete. {name} also contains a premature stop codon, which complicates the choice. The translated protein is: " + "-".join(protein)

    # Verify that Mutant 2 is the most disruptive
    # The protein from Mutant 2 is ['M', 'F', '_STOP_'], which is severely truncated.
    # The other mutations are missense or in-frame deletions, which are less likely
    # to cause a complete loss of function.
    
    # The LLM correctly identified Mutant 2 (Option D) as the one with the nonsense mutation.
    # Our analysis confirms this is the most likely mutation to cause loss of function.
    if llm_answer == "D":
        return "Correct"
    else:
        return f"Incorrect: The LLM chose {llm_answer}, but the correct answer is D (Mutant 2) because it is the only mutation that introduces a premature stop codon, leading to a non-functional truncated protein."

# Run the check
result = check_mutation_answer()
print(result)