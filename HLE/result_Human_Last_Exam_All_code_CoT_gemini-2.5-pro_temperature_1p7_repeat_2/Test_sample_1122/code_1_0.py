import io
from Bio.Blast import NCBIWWW, NCBIXML

def identify_protein_sequence(sequence: str):
    """
    Identifies a protein by performing a BLASTp search against the NCBI 'nr' database.

    Args:
        sequence (str): The amino acid sequence of the protein.
    """
    # Clean up the input sequence by removing spaces and newlines
    cleaned_sequence = "".join(sequence.split())
    
    print("Performing BLASTp search against the NCBI 'nr' database...")
    print("This may take a minute or two depending on server load.")
    
    try:
        # Use qblast to send the sequence to the NCBI BLAST web service
        # It returns a handle to the results in XML format
        result_handle = NCBIWWW.qblast(
            program="blastp",          # Protein-protein BLAST
            database="nr",             # Non-redundant protein database
            sequence=cleaned_sequence,
        )

        # To prevent potential issues with a closed handle, we read the result
        # into a string variable and then parse it.
        blast_xml_string = result_handle.read()
        result_handle.close()
        
        # Use StringIO to treat the XML string as a file for parsing
        blast_result_stream = io.StringIO(blast_xml_string)

        # Parse the XML results. We use read() because we only have one query.
        blast_record = NCBIXML.read(blast_result_stream)

        # Check if any alignments were found
        if blast_record.alignments:
            # The first alignment is the best match
            top_alignment = blast_record.alignments[0]
            top_hit = top_alignment.hsps[0] # High-scoring Segment Pair
            
            print("\n--- Search Complete ---")
            print(f"Top match found: {top_alignment.title}")
            print(f"Sequence Length: {top_alignment.length}")
            # E-value indicates the statistical significance of the match.
            # A value close to zero means the match is highly significant.
            print(f"E-value: {top_hit.expect:.2e}")
            
        else:
            print("\nNo significant matches found for the given sequence.")

    except Exception as e:
        print(f"\nAn error occurred during the BLAST search: {e}")
        print("Please check your internet connection or try again later, as the NCBI service might be busy.")


if __name__ == '__main__':
    # The amino acid sequence provided by the user
    amino_acid_sequence = """
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
    identify_protein_sequence(amino_acid_sequence)