import sys

# This script requires the 'biopython' library.
# If you don't have it, you can install it by running: pip install biopython
try:
    from Bio.Blast import NCBIWWW
    from Bio.Blast import NCBIXML
except ImportError:
    print("Error: The 'biopython' library is not installed.")
    print("Please install it by running the command: pip install biopython")
    sys.exit(1)

def identify_protein(sequence: str):
    """
    Identifies a protein sequence using NCBI BLAST.
    
    Args:
        sequence: The amino acid sequence as a string.
    """
    # Clean up the sequence by removing any whitespace or newlines
    cleaned_sequence = "".join(sequence.split())
    
    if not cleaned_sequence:
        print("The provided sequence is empty.")
        return

    print(f"Searching for protein with sequence starting: {cleaned_sequence[:30]}...")
    print("This may take a moment as it involves a live query to the NCBI database.")

    try:
        # Use qblast to query the NCBI 'blastp' program against the 'nr' database.
        result_handle = NCBIWWW.qblast("blastp", "nr", cleaned_sequence)
        
        # Parse the returned XML data.
        blast_record = NCBIXML.read(result_handle)
        
        # Check if any matches (alignments) were found.
        if blast_record.alignments:
            # The first alignment is the best match.
            top_hit = blast_record.alignments[0]
            top_hsp = top_hit.hsps[0] # Highest-scoring Pair
            
            # Extract relevant information
            protein_name = top_hit.title
            accession = top_hit.accession
            e_value = top_hsp.expect
            identity = (top_hsp.identities / top_hsp.align_length) * 100
            
            print("\n--- Protein Identification Complete ---")
            print(f"Best Match Found: {protein_name}")
            print(f"Accession ID: {accession}")
            print(f"E-value (Confidence): {e_value:.2e} (A value very close to 0 means a very high confidence match)")
            print(f"Identity: {identity:.2f}%")
            
        else:
            print("\nNo significant match was found for the provided sequence.")

    except Exception as e:
        print(f"\nAn error occurred during the search: {e}")
        print("Please check your internet connection and try again.")
    finally:
        if 'result_handle' in locals() and result_handle:
            result_handle.close()


if __name__ == "__main__":
    # The amino acid sequence to identify
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
    identify_protein(amino_acid_sequence)