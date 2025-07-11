import sys
# It's good practice to check if a third-party library is installed.
try:
    from Bio.Blast import NCBIWWW
    from Bio.Blast import NCBIXML
except ImportError:
    print("Biopython library not found.")
    print("Please install it using: pip install biopython")
    sys.exit()

def identify_protein(sequence: str):
    """
    Identifies a protein by performing a BLAST search against the NCBI database.
    
    Args:
        sequence (str): The amino acid sequence of the protein.
    """
    # Clean up the sequence string (remove spaces and newlines)
    cleaned_sequence = "".join(sequence.split())
    
    print("Performing BLASTP search against NCBI's 'nr' database...")
    print("This may take a moment. Please wait.")
    
    try:
        # Use qblast to submit the search.
        # program = "blastp" (protein-protein BLAST)
        # database = "nr" (non-redundant protein database)
        # sequence = our cleaned sequence
        result_handle = NCBIWWW.qblast("blastp", "nr", cleaned_sequence)
        
        # The result from NCBI is in XML format. We read it once.
        blast_record = NCBIXML.read(result_handle)
        
        # Check if any alignments were found
        if blast_record.alignments:
            # The first alignment is the best hit with the highest score
            first_hit = blast_record.alignments[0]
            # The hsp (high-scoring pair) contains statistics like the e-value
            hsp = first_hit.hsps[0]
            
            print("\n--- Top Match Found ---")
            print(f"Protein Name: {first_hit.title}")
            print(f"Accession ID: {first_hit.accession}")
            print(f"Sequence Length: {first_hit.length} aa")
            print(f"E-value (Expect value): {hsp.expect}")
            print(f"Score: {hsp.score}")
            print(f"Identities: {hsp.identities}/{hsp.align_length} ({int(hsp.identities/hsp.align_length*100)}%)")
            print("-----------------------")
        else:
            print("\nNo significant matches were found for the provided sequence.")
            
    except Exception as e:
        print(f"\nAn error occurred during the BLAST search: {e}")
        print("Please check your internet connection and try again.")
    finally:
        # It's important to close the handle
        if 'result_handle' in locals() and result_handle:
            result_handle.close()

if __name__ == '__main__':
    # The amino acid sequence provided by the user
    input_sequence = """
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
    identify_protein(input_sequence)