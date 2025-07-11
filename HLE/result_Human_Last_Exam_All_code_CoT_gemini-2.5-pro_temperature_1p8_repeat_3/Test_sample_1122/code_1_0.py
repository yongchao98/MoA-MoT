import io
from urllib.error import URLError

try:
    from Bio.Blast import NCBIWWW, NCBIXML
except ImportError:
    print("Error: Biopython is not installed.")
    print("Please install it first by running: pip install biopython")
    exit()

def identify_protein_from_sequence():
    """
    Identifies a protein by performing a BLAST search against the NCBI database.
    """
    # The amino acid sequence to be identified.
    sequence_string = (
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

    print("Connecting to NCBI and performing a BLAST search...")
    print("This may take a few moments depending on server load.")

    try:
        # Submit the BLAST search via NCBI's web API.
        # 'blastp' is for protein-protein search. 'nr' is the non-redundant protein db.
        result_handle = NCBIWWW.qblast("blastp", "nr", sequence_string)

        # The result comes in XML format. We read and parse it.
        # We use io.StringIO to treat the XML string result as a file.
        blast_xml_string = result_handle.read()
        result_handle.close()
        blast_in = io.StringIO(blast_xml_string)

        # Use NCBIXML.parse to read one or more BLAST records.
        # Since we sent one query, we expect one record.
        blast_records = list(NCBIXML.parse(blast_in))

        if not blast_records:
             print("\nNo BLAST records were parsed. The server might be busy or the result was empty.")
             return
             
        blast_record = blast_records[0]

        # Check if the search returned any matching alignments.
        if blast_record.alignments:
            # The first alignment is the best hit.
            top_hit = blast_record.alignments[0]
            top_hsp = top_hit.hsps[0] # High-scoring Segment Pair

            print("\n--- Search Complete: Top Match Found ---")
            print(f"Description: {top_hit.title}")
            print(f"Sequence Length: {top_hit.length}")
            print(f"E-value (Expect): {top_hsp.expect}")
            print(f"Identity: {top_hsp.identities}/{top_hsp.align_length} ({100 * top_hsp.identities / top_hsp.align_length:.2f}%)")

            # Extract a concise answer for the final output
            # e.g. from ">ref|NP_000658.1| alcohol dehydrogenase 1A (class I), alpha polypeptide [Homo sapiens]"
            title_parts = top_hit.title.split('[')
            protein_name = title_parts[0].split('|')[-1].strip()
            species = title_parts[1][:-1].strip() if len(title_parts) > 1 else "Unknown species"
            print(f"\nConclusion: The sequence represents '{protein_name}' from the organism '{species}'.")
        else:
            print("\nNo significant matches were found for the provided sequence.")

    except URLError:
        print("\nNetwork Error: Could not connect to the NCBI BLAST server.")
        print("Please check your internet connection.")
    except Exception as e:
        print(f"\nAn unexpected error occurred: {e}")

if __name__ == "__main__":
    identify_protein_from_sequence()
