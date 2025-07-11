import io
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

def identify_protein_sequence():
    """
    Identifies a protein by performing a BLASTp search against the NCBI 'nr' database.
    """
    # The amino acid sequence to be identified
    sequence_data = (
        "MAEQVALSRTQVCGILREELFQGDAFHQSDTHIFIIMGASGDLAKKKIYP"
        "TIWWLFRDGLLPENTFIVGYARSRLTVADIRKQSEPFFKATPEEKLKLED"
        "FFARNSYVAGQYDDAASYQRLNSHMNALHLGSQANRLFYLALPPTVYEAV"
        "TKNIHESCMQIGWNRIIVEKPFGRDLQSSDRLSNHISSLFREDQIYRIDH"
        "YLGKEMVQNLMVLRFANRIFGPIWNRDNIACVILTFKEPFGTEGRGGYFD"
        "EFGIIRDVMQNHLLQMLCLVAMEKPASTNSDDVRDEKVKVLKCISEVQAN"
        "NVVLGQYVGNPDGEGEATKGYLDDPTVPRGSTATTFAAVVLYVENERWDG"
        "VPFILRCGKALNERKAEVRLQFHDVAGDIFHQQCKRNELVIRVQPNEAVY"
        "TKMMTKKPGMFFNPEESELDLTYGNRYKNVKLPDAYERLILDVFCGSQMH"
        "FVRSDELREAWRIFTPLLHQIELEKPKPIPYIYGSRGPTEADELMKRVGF"
        "QYEGTYKWVNPHKL"
    )

    print("Submitting sequence to NCBI BLAST. This may take a moment...")

    try:
        # Use qblast to query the NCBI BLAST server
        # program='blastp' for protein-protein search
        # database='nr' for the non-redundant protein database
        result_handle = NCBIWWW.qblast("blastp", "nr", sequence_data)

        # The result from NCBI is in XML format. We need to parse it.
        # We read the result into memory before parsing.
        blast_xml_result = result_handle.read()
        result_handle.close()
        
        # Use StringIO to treat the XML string as a file
        blast_in = io.StringIO(blast_xml_result)
        blast_record = NCBIXML.read(blast_in)

        # Check if any alignments were found
        if blast_record.alignments:
            # The top alignment is the most likely match
            top_alignment = blast_record.alignments[0]
            print("\n--- Top Match Found ---")
            print(f"Description: {top_alignment.title}")
            
            # The hsp (high-scoring pair) object contains statistics
            top_hsp = top_alignment.hsps[0]
            print(f"Accession ID: {top_alignment.accession}")
            print(f"Sequence Length: {top_alignment.length}")
            print(f"E-value: {top_hsp.expect}")
            print(f"Score: {top_hsp.score}")
            print(f"Identities: {top_hsp.identities}/{top_hsp.align_length} ({100 * top_hsp.identities / top_hsp.align_length:.2f}%)")

        else:
            print("No significant match was found for the given sequence.")

    except Exception as e:
        print(f"An error occurred during the BLAST search: {e}")
        print("Please ensure you have an active internet connection and that the biopython library is installed.")

if __name__ == "__main__":
    identify_protein_sequence()