# First, ensure you have biopython installed:
# pip install biopython

import sys
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

def identify_protein(sequence):
    """
    Identifies a protein sequence using NCBI BLAST.
    """
    print("Performing BLAST search online. This may take a few moments...")
    
    try:
        # Use NCBI's online BLAST service (qblast)
        # 'blastp' is for protein-protein comparison
        # 'nr' is the non-redundant protein database
        result_handle = NCBIWWW.qblast("blastp", "nr", sequence)
        
        # The result from qblast is in XML format. We need to parse it.
        # We use read() because we performed a single query.
        blast_record = NCBIXML.read(result_handle)
        
        # Check if any results were returned
        if not blast_record.alignments:
            print("No significant matches found for the provided sequence.")
            return

        # The first alignment is the best hit
        top_alignment = blast_record.alignments[0]
        top_hsp = top_alignment.hsps[0] # High-scoring Segment Pair

        # Extracting relevant information
        protein_title = top_alignment.title
        sequence_length = top_alignment.length
        e_value = top_hsp.expect
        identities = top_hsp.identities
        alignment_length = top_hsp.align_length
        identity_percentage = (identities / alignment_length) * 100

        print("\n--- Top Match Found ---")
        print(f"Protein Name: {protein_title}")
        print(f"Sequence Length: {sequence_length}")
        print(f"E-value: {e_value}")
        
        # The user requested to output each number.
        # Here we format it like an "equation" for identity calculation.
        print("Identity Calculation:")
        print(f"{identities} (Identical Residues) / {alignment_length} (Alignment Length) = {identity_percentage:.2f}% Identity")

    except Exception as e:
        print(f"An error occurred: {e}", file=sys.stderr)
        print("Please check your internet connection and ensure the NCBI services are available.", file=sys.stderr)
    finally:
        if 'result_handle' in locals() and result_handle:
            result_handle.close()

if __name__ == "__main__":
    # The amino acid sequence to identify
    amino_acid_sequence = (
        "MAEQVALSRTQVCGILREELFQGDAFHQSDTHIFIIMGASGDLAKKKIYP"
        "TIWWLFRDGLLPENTFIVGYARSRLTVADIRKQSEPFFKATPEEKLKLED"
        "FFARNSYVAGQYDDAASYQRLNSHMNALHLGSQANRLFYLALPPTVYEAV"
        "TKNIHESCMQIGWNRIIVEKPFGRDLQSSDRLSNHISSLFREDQIYRIDH"
        "YLGKEMVQNLMVLRFANRIFGPIWNRDNIACVILTFKEPFGTEGRGGYFD"
        "EFGIIRDVMQNHLLQMLCLVAMEKPASTNSDDVRDEKVKVLKCISEVQAN"
        "NVVLGQYVGNPDGEGEATKGYLDDPTVPRGSTTATFAAVVLYVENERWDG"
        "VPFILRCGKALNERKAEVRLQFHDVAGDIFHQQCKRNELVIRVQPNEAVY"
        "TKMMTKKPGMFFNPEESELDLTYGNRYKNVKLPDAYERLILDVFCGSQMH"
        "FVRSDELREAWRIFTPLLHQIELEKPKPIPIYIGSRGPTEADELMKRVGF"
        "QYEGTYKWVNPHKL"
    )
    
    # Remove any spaces or newlines from the sequence string
    cleaned_sequence = "".join(amino_acid_sequence.split())
    
    identify_protein(cleaned_sequence)