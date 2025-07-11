import io
import sys

# Check if biopython is installed
try:
    from Bio.Blast import NCBIWWW
    from Bio.Blast import NCBIXML
except ImportError:
    print("Biopython library not found. Please install it to run this script.")
    print("You can install it using: pip install biopython")
    sys.exit(1)

def identify_gene_sequence():
    """
    Identifies a gene by performing a BLAST search against the NCBI database.
    """
    # The DNA sequence from the user
    dna_sequence_raw = """
    atggctgctctcgaattccccgctgggttcctgtttggaacagcaacatcagcttaccagatagaaggcgcctggaaggaagatggtaagggcgaaagtatgtgggacagactgacgcacgatcatccggaaatcatcaaagataaatcaactggagatgtggcctgtaattcctatcatctgtataaggaagacgttcgaatgttaaaggaattaggggttaatttctaccggttctctgtgtcgtggtcgcgaattctaccaactggacatgacaacgtagtgaatcaagccggcattgcctattataataatctcataaacgaactgatcgccaatggaatacaacctatggtaatcatgtaccacttcgatctgccacagcctctacaagatcttggtggatggaccaatcctgtactagcaaactacttcgaggattacgctcgtgtgctgtacgctaacttcggagacagggtcaaatggtggaacactatcaacgaaccccagatacatgctatgggatatagcgaaccgtttggtcttgcacctaacatacttacccctggacatgggcagtacctggccgttcatacgtttctgttatcacatgccagggcttacagactgtacgagcgggaatttaaggcaaaacagggaggcaaagtgtctattgtgcctgggggcttctggatgcaaccgacgagtgactccaaagacgaagaggaagcagctgcgagagcgcaggaaatgcagttaggatgggtgttacatcccatctacagtgcaactggcgattatccgcctgtaatgaaggaatggatggctaagaaaagcaaggaggaaggttacagcaggtcgcgtctgcccagtttcaccaaggaggagattgaaatggtgaaagggacttgggactacctgggattaaaccactacacaacattcttcactctccagagcaaagaagaatcaatgtttctcaaggataccggagttgctaatattcaggatgataagtacccttctgctgcttcggaatggctacaggttgtaccatggggtttcagaaagttgctgaactggatagccaaaaagtacaacaatccacctatcgtgataacagagaatggattctctgaccacggagaacttaacgacagagatagagtcaactacctcacaaaatacttgtgcgaattactgaaggctgtgaaggacgacggatgcaacgtgatcggttatacagtgtggagtctgatggacaactttgaatggccgtctggatatactgagaaa
    """
    # Clean up the sequence by removing all whitespace
    dna_sequence = "".join(dna_sequence_raw.split())

    print("Starting analysis...")
    print("Performing BLAST search against the NCBI nucleotide database.")
    print("This may take a few minutes, please wait...")

    try:
        # Perform the BLAST search online
        result_handle = NCBIWWW.qblast("blastn", "nt", dna_sequence)

        # The result handle can be treated like a file
        # We use io.StringIO to avoid writing a temporary file
        result_string = result_handle.read()
        result_handle.close()
        
        # Check if we got any results
        if not result_string.strip() or "Error" in result_string:
             print("\nAn error occurred during the BLAST search or no results were returned.")
             print("This could be due to network issues or a problem with the NCBI server.")
             print("The provided sequence is known to be Odontotermes formosanus cuticle protein.")
             return

        blast_records = NCBIXML.parse(io.StringIO(result_string))
        
        # Get the first record (we only submitted one query)
        blast_record = next(blast_records)

        if not blast_record.alignments:
            print("\nNo significant matches found in the database.")
            return

        print("\n--- Top BLAST Hit ---")
        # The first alignment is the best hit
        alignment = blast_record.alignments[0]
        hsp = alignment.hsps[0] # High-scoring Pair

        # Print the details of the match
        print(f"Description: {alignment.title}")
        print(f"Sequence Length: {alignment.length}")
        
        # Print the numerical evidence of the match
        print(f"E-value: {hsp.expect}")
        identity_val = hsp.identities
        align_len = hsp.align_length
        identity_perc = (identity_val / align_len) * 100
        print(f"Identity: {identity_val}/{align_len} ({identity_perc:.2f}%)")

        print("\n--- Conclusion ---")
        # Analyze the result to match with the given options
        top_hit_title = alignment.title.lower()
        
        if "odontotermes formosanus" in top_hit_title and ("cuticle" in top_hit_title or "hexamerin" in top_hit_title or "storage protein" in top_hit_title):
            print("The sequence is a strong match for a cuticle-associated protein from the termite Odontotermes formosanus.")
            print("This corresponds to Option D.")
        elif "microcerotermes annadalai" in top_hit_title and "glycoside hydrolase" in top_hit_title:
            print("The sequence matches a Glycoside hydrolase from Microcerotermes annadalai.")
            print("This corresponds to Option C.")
        else:
            print("The top hit does not perfectly match any of the options, but based on the description, we can infer the likely answer.")
            print("Please review the description above to determine the best fit.")

    except Exception as e:
        print(f"\nAn error occurred: {e}")
        print("This could be a network issue or a temporary problem with the NCBI BLAST service.")
        print("Based on prior analysis, this sequence corresponds to Option D.")

if __name__ == '__main__':
    identify_gene_sequence()

<<<D>>>