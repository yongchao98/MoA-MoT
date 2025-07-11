# First, ensure you have biopython installed:
# pip install biopython

from Bio.Blast import NCBIWWW, NCBIXML

def identify_protein_sequence(sequence: str):
    """
    Identifies a protein sequence using NCBI's online BLAST service.

    Args:
        sequence: A string representing the amino acid sequence.
    """
    print("Performing BLAST search against the NCBI database...")
    print("This may take a minute or two...")

    try:
        # Perform the BLAST search using qblast
        # program='blastp' for protein-protein search
        # database='nr' for the non-redundant protein database
        result_handle = NCBIWWW.qblast("blastp", "nr", sequence)

        # Parse the XML result from BLAST
        blast_record = NCBIXML.read(result_handle)

        # Check if any alignments were found
        if not blast_record.alignments:
            print("No matching protein found for the given sequence.")
            return

        # The first alignment is the best match
        top_alignment = blast_record.alignments[0]
        top_hit = top_alignment.hsps[0] # High-scoring Segment Pair

        print("\n--- Top Match Found ---")
        # The title of the alignment contains the protein name and species
        print(f"Description: {top_alignment.title}")
        print(f"Length: {top_alignment.length}")
        # The E-value (Expect value) represents the number of hits one can "expect"
        # to see by chance. A very low E-value indicates a significant match.
        print(f"E-value: {top_hit.expect}")
        print(f"Identity: {top_hit.identities}/{top_hit.align_length} ({100 * top_hit.identities / top_hit.align_length:.2f}%)")

        # Extract just the protein name for the final answer
        # The title is often formatted like 'gi|...|ref|NP_...| protein name [Species]'
        # We can try to clean it up for a concise answer.
        protein_info = top_alignment.title.split('|')[-1].strip()
        print(f"\nThe sequence most likely represents: {protein_info}")


    except Exception as e:
        print(f"An error occurred during the BLAST search: {e}")
    finally:
        if 'result_handle' in locals() and result_handle:
            result_handle.close()

if __name__ == "__main__":
    # The amino acid sequence provided by the user
    amino_acid_sequence = (
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

    # Remove spaces and newlines to create a single sequence string
    cleaned_sequence = "".join(amino_acid_sequence.split())

    identify_protein_sequence(cleaned_sequence)
