import io
from Bio.Blast import NCBIWWW, NCBIXML

def identify_protein_sequence():
    """
    Identifies a protein by performing a BLASTp search with its amino acid sequence.
    """
    # The amino acid sequence provided by the user
    sequence_data = """
    MAEQVALSRT QVCGILREEL FQGDAFHQSD THIFIIMGAS GDLAKKKIYP 
    TIWWLFRDGL LPENTFIVGY ARSRLTVADI RKQSEPFFKA TPEEKLKLED 
    FFARNSYVAG QYDDAASYQR LNSHMNALHL GSQANRLFYL ALPPTVYEAV 
    TKNIHESCMQ IGWNRIIVEK PFGRDLQSSD RLSNHISSLF REDQIYRIDH 
    YLGKEMVQNL MVLRFANRIF GPIWNRDNIA CVILTFKEPF GTEGRGGYFD 
    EFGIIRDVMQ NHLLQMLCLV AMEKPASTNS DDVRDEKVKV LKCISEVQAN 
    NVVLGQYVGN PDGEGEATKG YLDDPTVPRG STTATFAAVV LYVENERWDG 
    VPFILRCGKA LNERKAEVRL QFHDVAGDIF HQQCKRNELV IRVQPNEAVY 
    TKMMTKKPGM FFNPEESELD LTYGNRYKNV KLPDAYERLI LDVFCGSQMH 
    FVRSDELREA WRIFTPLLHQ IELEKPKPIP YIYGSRGPTE ADELMKRVGF 
    QYEGTYKWVN PHKL
    """
    
    # Clean up the sequence by removing whitespace and newlines
    sequence = "".join(sequence_data.strip().split())
    
    print("Performing a BLASTp search against the NCBI 'nr' database.")
    print("This may take a minute or two, please wait...")

    try:
        # Perform the BLAST search over the internet
        result_handle = NCBIWWW.qblast("blastp", "nr", sequence)
        
        # The result comes back as XML. We read it into a string to parse it.
        # This is done to prevent issues with reading a handle multiple times.
        blast_xml_string = result_handle.read()
        result_handle.close()
        
        # Use io.StringIO to treat the XML string as a file for the parser
        blast_records_handle = io.StringIO(blast_xml_string)
        
        # The NCBIXML.parse() function returns an iterator of blast records.
        # For a single query, there will be only one record.
        blast_record = next(NCBIXML.parse(blast_records_handle))

        # Check if any alignments were found
        if not blast_record.alignments:
            print("\nNo significant matches found for the given sequence.")
        else:
            # The first alignment is the best hit
            top_alignment = blast_record.alignments[0]
            # The first HSP (High-scoring Segment Pair) contains score details
            top_hsp = top_alignment.hsps[0]
            
            print("\n--- Top Match Found ---")
            # The 'title' contains the full description of the protein match
            print(f"Description: {top_alignment.title}")
            print(f"Accession ID: {top_alignment.accession}")
            # The E-value (Expect value) represents the number of hits one can "expect" 
            # to see by chance when searching a database of a particular size.
            # A lower E-value indicates a more significant match.
            print(f"E-value: {top_hsp.expect}")
            print(f"Score: {top_hsp.score}")
            print(f"Identities: {top_hsp.identities}/{top_hsp.align_length} ({100.0 * top_hsp.identities / top_hsp.align_length:.2f}%)")

    except Exception as e:
        print(f"\nAn error occurred: {e}")
        print("This could be due to a network issue or a problem with the NCBI servers.")
        print("Please ensure you have an active internet connection and have the 'biopython' library installed (`pip install biopython`).")

if __name__ == "__main__":
    identify_protein_sequence()