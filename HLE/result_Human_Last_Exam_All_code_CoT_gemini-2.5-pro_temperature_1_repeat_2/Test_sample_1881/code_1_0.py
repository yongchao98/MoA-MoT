import sys
# First, ensure biopython is installed
try:
    from Bio.Blast import NCBIWWW
    from Bio.Blast import NCBIXML
except ImportError:
    print("Biopython library not found.")
    print("Please install it using: pip install biopython")
    sys.exit(1)

def identify_gene_sequence():
    """
    Performs a BLAST search on a given DNA sequence to identify its function
    and origin.
    """
    # The DNA sequence provided in the problem
    dna_sequence = """
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
    # Clean up the sequence by removing newlines and spaces
    cleaned_sequence = "".join(dna_sequence.split())

    print("Performing BLASTX search against the NCBI 'nr' protein database.")
    print("This may take a minute or two, please wait...")

    result_handle = None
    try:
        # Use blastx to compare the translated DNA sequence against a protein database
        result_handle = NCBIWWW.qblast("blastx", "nr", cleaned_sequence)

        # Parse the XML result from BLAST
        blast_record = NCBIXML.read(result_handle)

        # Print the results header
        print("\n--- Top BLAST Hits ---")
        print(f"{'Description':<80} {'E-value':<10}")
        print("-" * 90)

        # Check if there are any hits
        if not blast_record.alignments:
            print("No significant matches found.")
            return

        # Display the top 5 hits
        limit = 5
        for i, alignment in enumerate(blast_record.alignments):
            if i >= limit:
                break
            # The E-value (expect value) is for the first and best high-scoring pair (HSP)
            e_value = alignment.hsps[0].expect
            print(f"{alignment.title[:78]:<80} {f'{e_value:.2e}':<10}")
            
        print("\n--- Conclusion ---")
        top_hit_title = blast_record.alignments[0].title
        if "glycoside hydrolase" in top_hit_title.lower() and \
           ("microcerotermes" in top_hit_title.lower() or "odontotermes" in top_hit_title.lower()):
            print("The results strongly indicate the protein is a Glycoside Hydrolase from a termite species.")
            print("Comparing this with the given options, 'C' is the best match.")
        else:
            print("The top hit does not conclusively match any of the options, but analysis of the top hit is the correct methodology.")


    except Exception as e:
        print(f"\nAn error occurred during the BLAST search: {e}")
        print("This may be due to a network issue or temporary unavailability of the NCBI service.")
    finally:
        # It's good practice to close the handle
        if result_handle:
            result_handle.close()

if __name__ == "__main__":
    identify_gene_sequence()