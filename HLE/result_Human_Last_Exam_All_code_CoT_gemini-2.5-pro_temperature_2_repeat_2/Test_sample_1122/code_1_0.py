import sys
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

def identify_protein(sequence):
    """
    Identifies a protein by performing a BLAST search with its amino acid sequence.

    Args:
        sequence (str): The amino acid sequence of the protein.
    """
    # Clean up the sequence by removing any spaces or newline characters
    clean_sequence = "".join(sequence.split())

    try:
        print("Performing a BLAST search against the NCBI database.")
        print("This may take a moment...")
        # Use NCBI's online BLAST tool (blastp for protein-protein search)
        result_handle = NCBIWWW.qblast("blastp", "nr", clean_sequence)

        # Parse the XML output from the BLAST search
        blast_record = NCBIXML.read(result_handle)

        # Check if any alignments (matches) were found
        if blast_record.alignments:
            # The first alignment is the best hit
            top_hit = blast_record.alignments[0]
            print("\n--- Top Match Found ---")
            # The 'title' contains the description of the matching protein
            print(f"Description: {top_hit.title}")
            # The 'length' is the number of amino acids in the matching sequence
            print(f"Length: {top_hit.length}")

            # Also show some statistics of the best alignment segment (HSP)
            if top_hit.hsps:
                top_hsp = top_hit.hsps[0]
                print(f"E-value: {top_hsp.expect}")
                print(f"Score: {top_hsp.score}")
                print(f"Identity: {top_hsp.identities}/{top_hsp.align_length} ({(top_hsp.identities/top_hsp.align_length)*100:.2f}%)")
        else:
            print("No significant match was found for the provided sequence.")

    except Exception as e:
        print(f"An error occurred: {e}", file=sys.stderr)
        print("Please check your internet connection and ensure the sequence is valid.", file=sys.stderr)
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