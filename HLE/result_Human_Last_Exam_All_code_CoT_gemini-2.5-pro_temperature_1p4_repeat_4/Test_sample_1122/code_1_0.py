import sys

# A helper function to check for and guide the user to install a missing library.
def install_and_import(package):
    import importlib
    try:
        importlib.import_module(package)
    except ImportError:
        print(f"The '{package}' library is required to run this script.")
        print(f"Please install it by running: pip install {package}")
        sys.exit(1)
    finally:
        globals()[package] = importlib.import_module(package)

# Ensure biopython is available
install_and_import('Bio')
from Bio.Blast import NCBIWWW, NCBIXML


def identify_protein(sequence: str):
    """
    Identifies a protein sequence using NCBI's BLAST service.
    
    Args:
        sequence: The amino acid sequence of the protein.
    """
    print("Submitting sequence to NCBI BLAST for identification...")
    print("This may take a minute or two, please wait.")
    
    try:
        # Use qblast to perform an online BLAST search.
        # "blastp" is for comparing a protein sequence against a protein database.
        # "nr" is the non-redundant protein database, which is comprehensive.
        result_handle = NCBIWWW.qblast("blastp", "nr", sequence)
        
        # The result from NCBI is in XML format. We use NCBIXML to parse it.
        blast_record = NCBIXML.read(result_handle)
        
        # Check if any alignments (matches) were found.
        if blast_record.alignments:
            # The first alignment in the record is the best match.
            top_hit = blast_record.alignments[0]
            print("\n--- Top Match Found ---")
            print(f"Protein: {top_hit.title}")
            print(f"Length: {top_hit.length} amino acids")
            # The E-value (Expect value) represents the number of hits one can "expect" 
            # to see by chance. A smaller E-value indicates a more significant match.
            print(f"E-value: {top_hit.hsps[0].expect}")
        else:
            print("\nNo significant match was found for the provided sequence.")
            
    except Exception as e:
        print(f"\nAn error occurred during the BLAST search: {e}")
        print("This could be a temporary network issue or a problem with the NCBI server.")
        print("Please check your internet connection and try again later.")
        
    finally:
        # It's good practice to close the handle.
        if 'result_handle' in locals() and result_handle:
            result_handle.close()

# The amino acid sequence from the user.
protein_sequence = "MAEQVALSRTQVCGILREELFQGDAFHQSDTHIFIIMGASGDLAKKKIYPTIWWLFRDGLLPENTFIVGYARSRLTVADIRKQSEPFFKATPEEKLKLEDFFARNSYVAGQYDDAASYQRLNSHMNALHLGSQANRLFYLALPPTVYEAVTKNIHESCMQIGWNRIIVEKPFGRDLQSSDRLSNHISSLFREDQIYRIDHYLGKEMVQNLMVLRFANRIFGPIWNRDNIACVILTFKEPFGTEGRGGYFDEFGIIRDVMQNHLLQMLCLVAMEKPASTNSDDVRDEKVKVLKCISEVQANNVVLGQYVGNPDGEGEATKGYLDDPTVPRGSTATTFAAVVLYVENERWDGVPFILRCGKALNERKAEVRLQFHDVAGDIFHQQCKRNELVIRVQPNEAVYTMMTKKPGMFFNPEESELDLTYGNRYKNVKLPDAYERLILDVFCGSQMHFVRSDELREAWRIFTPLLHQIELEKPKPIPIYGSRGPTEADELMKRVGFQYEGTYKWVNPHKL"

identify_protein(protein_sequence)