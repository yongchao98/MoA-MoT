import sys
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

def identify_protein_sequence(sequence):
    """
    Performs a BLASTP search for a given protein sequence and prints the top hit.
    """
    # Clean up the sequence string by removing spaces and newlines
    cleaned_sequence = "".join(sequence.split())
    
    print("Performing BLAST search against the 'nr' database. This may take a few moments...")
    
    try:
        # Use the NCBIWWW.qblast function for the search
        # 'blastp' is for protein-protein comparison
        # 'nr' is the non-redundant protein database
        result_handle = NCBIWWW.qblast("blastp", "nr", cleaned_sequence)
        
        # The result from NCBI is in XML format, which we need to parse
        blast_record = NCBIXML.read(result_handle)
        
        # Check if any results were found
        if not blast_record.alignments:
            print("No matching protein found for the provided sequence.")
            return "No protein found"
            
        # The first alignment is the best match
        top_alignment = blast_record.alignments[0]
        top_hit_title = top_alignment.title
        
        # The HSP (High-scoring Segment Pair) contains statistics about the match
        top_hsp = top_alignment.hsps[0]
        
        print("\n--- Top Match Found ---")
        print(f"Description: {top_hit_title}")
        print(f"Sequence Length: {top_alignment.length}")
        print(f"E-value (Expect value): {top_hsp.expect}")
        print(f"Score: {top_hsp.score}")
        
        # The final answer is the description of the top hit
        final_answer = top_hit_title

    except Exception as e:
        print(f"An error occurred during the BLAST search: {e}", file=sys.stderr)
        print("This could be due to a network issue or a problem with the NCBI service.", file=sys.stderr)
        print("Please ensure you have an internet connection and the biopython library is installed ('pip install biopython').", file=sys.stderr)
        final_answer = "Error during search"
        
    finally:
        # Always close the handle
        if 'result_handle' in locals():
            result_handle.close()
            
    return final_answer

if __name__ == '__main__':
    # The amino acid sequence provided by the user
    user_sequence = (
        "MAEQVALSRT QVCGILREEL FQGDAFHQSD THIFIIMGAS GDLAKKKIYP "
        "TIWWLFRDGL LPENTFIVGY ARSRLTVADI RKQSEPFFKA TPEEKLKLED "
        "FFARNSYVAG QYDDAASYQR LNSHMNALHL GSQANRLFYL ALPPTVYEAV "
        "TKNIHESCMQ IGWNRIIVEK PFGRDLQSSD RLSNHISSLF REDQIYRIDH "
        "YLGKEMVQNL MVLRFANRIF GPIWNRDNIA CVILTFKEPF GTEGRGGYFD "
        "EFGIIRDVMQ NHLLQMLCLV AMEKPASTNS DDVRDEKVKV LKCISEVQAN "
        "NVVLGQYVGN PDGEGEATKG YLDDPTVPRG STTATFAAVV LYVENERWDG "
        "VPFILRCGKA LNERKAEVRL QFHDVAGDIF HQQCKRNELV IRVQPNEAVY "
        "TKMMTKKPGM FFNPEESELD LTYGNRYKNV KLPDAYERLI LDVFCGSQMH "
        "FVRSDELREA WRIFTPLLHQ IELEKPKPIP YIYGSRGPTE ADELMKRVGF "
        "QYEGTYKWVN PHKL"
    )
    
    # Run the identification function
    result = identify_protein_sequence(user_sequence)
