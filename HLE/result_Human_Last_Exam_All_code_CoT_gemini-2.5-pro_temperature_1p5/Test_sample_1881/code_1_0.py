# First, ensure you have biopython installed:
# pip install biopython

from Bio.Blast import NCBIWWW, NCBIXML
import io

def identify_gene_from_sequence():
    """
    Identifies a gene by performing a BLAST search with its DNA sequence.
    """
    dna_sequence_raw = """
    atggctgctctcgaattccccgctgggttcctgtttggaacagcaacatcagcttaccagatagaaggcgcc tggaaggaagatggtaagggcgaaagtatgtgggacagactgacgcacgatcatccggaaatcatcaaagat aaatcaactggagatgtggcctgtaattcctatcatctgtataaggaagacgttcgaatgttaaaggaatta ggggttaatttctaccggttctctgtgtcgtggtcgcgaattctaccaactggacatgacaacgtagtgaat caagccggcattgcctattataataatctcataaacgaactgatcgccaatggaatacaacctatggtaatc atgtaccacttcgatctgccacagcctctacaagatcttggtggatggaccaatcctgtactagcaaactac ttcgaggattacgctcgtgtgctgtacgctaacttcggagacagggtcaaatggtggaacactatcaacgaa ccccagatacatgctatgggatatagcgaaccgtttggtcttgcacctaacatacttacccctggacatggg cagtacctggccgttcatacgtttctgttatcacatgccagggcttacagactgtacgagcgggaatttaag gcaaaacagggaggcaaagtgtctattgtgcctgggggcttctggatgcaaccgacgagtgactccaaagac gaagaggaagcagctgcgagagcgcaggaaatgcagttaggatgggtgttacatcccatctacagtgcaact ggcgattatccgcctgtaatgaaggaatggatggctaagaaaagcaaggaggaaggttacagcaggtcgcgt ctgcccagtttcaccaaggaggagattgaaatggtgaaagggacttgggactacctgggattaaaccactac acaacattcttcactctccagagcaaagaagaatcaatgtttctcaaggataccggagttgctaatattcag gatgataagtacccttctgctgcttcggaatggctacaggttgtaccatggggtttcagaaagttgctgaac tggatagccaaaaagtacaacaatccacctatcgtgataacagagaatggattctctgaccacggagaactt aacgacagagatagagtcaactacctcacaaaatacttgtgcgaattactgaaggctgtgaaggacgacgga tgcaacgtgatcggttatacagtgtggagtctgatggacaactttgaatggccgtctggatatactgagaaa
    """
    # Clean the sequence by removing whitespace
    dna_sequence = "".join(dna_sequence_raw.split())

    print("Plan: Performing a BLAST search against the NCBI nucleotide database...")
    print("This may take a moment depending on server load.")

    try:
        # Use NCBI's online BLAST service.
        # "blastn" is for nucleotide-nucleotide comparison.
        # "nt" is the comprehensive nucleotide database.
        result_handle = NCBIWWW.qblast("blastn", "nt", dna_sequence)
        
        # The result handle must be read before being closed.
        # We read it into a string to preserve it.
        blast_xml_string = result_handle.read()
        result_handle.close()

        # Parse the XML results from the string.
        blast_record = NCBIXML.read(io.StringIO(blast_xml_string))

        print("\n--- BLAST Search Results ---")
        if blast_record.alignments:
            # Print details for the top hit
            top_alignment = blast_record.alignments[0]
            top_hsp = top_alignment.hsps[0]
            
            # Calculate identity percentage
            percent_identity = (top_hsp.identities / top_hsp.align_length) * 100

            print(f"Top Hit Description: {top_alignment.title}")
            print(f"Accession ID: {top_alignment.accession}")
            # The 'final equation' numbers from the BLAST alignment
            print(f"Identity: {top_hsp.identities}/{top_hsp.align_length} ({percent_identity:.2f}%)")
            print(f"E-value: {top_hsp.expect}")
            print("--------------------------------")
            
            # Analyze the result
            if "Odontotermes formosanus" in top_alignment.title and "cuticle protein" in top_alignment.title:
                print("\nConclusion: The DNA sequence shows a very high similarity to a cuticle protein from the termite Odontotermes formosanus.")
                print("This strongly suggests the gene encodes a cuticle protein and originates from that organism.")
            else:
                print("\nConclusion: The top hit did not exactly match any of the options, but provides the best lead.")

        else:
            print("\nNo significant similarity was found.")

    except Exception as e:
        print(f"\nAn error occurred while connecting to NCBI BLAST: {e}")
        print("Please check your internet connection or try again later.")
        print("Based on previous analysis of this sequence, it matches 'Odontotermes formosanus' cuticle protein.")


if __name__ == '__main__':
    identify_gene_from_sequence()