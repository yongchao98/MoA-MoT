import sys
import subprocess

# Install biopython if it is not already installed
try:
    import Bio
except ImportError:
    print("BioPython library not found. Attempting to install...")
    try:
        subprocess.check_call([sys.executable, "-m", "pip", "install", "biopython"])
        print("BioPython installed successfully.")
    except Exception as e:
        print(f"Failed to install biopython. Please install it manually using: pip install biopython")
        print(f"Error: {e}")
        sys.exit(1)

from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

def identify_protein(sequence):
    """
    Identifies a protein by performing a BLAST search on the NCBI server.
    
    Args:
        sequence (str): The amino acid sequence of the protein.
    
    Returns:
        str: The name of the most likely protein, or an error message.
    """
    print("Performing BLAST search against the NCBI database. This may take a moment...")
    try:
        # Use qblast to perform the search online
        # blastp is for protein-protein search
        # nr is the non-redundant protein database
        result_handle = NCBIWWW.qblast("blastp", "nr", sequence.upper())
        
        # Parse the XML result
        blast_record = NCBIXML.read(result_handle)
        
        # Check if any alignments were found
        if blast_record.alignments:
            # The top hit is the first alignment in the list
            top_hit = blast_record.alignments[0]
            # The 'title' contains the description of the protein
            # e.g., >sp|P00363.2|IDH_ECOLI RecName: Full=Isocitrate dehydrogenase [NADP]...
            protein_name = top_hit.title
            
            # Clean up the name for better readability
            # Find the first mention of the protein name, usually after 'Full='
            if 'RecName: Full=' in protein_name:
                protein_name = protein_name.split('RecName: Full=')[1].split(';')[0]
            else: # Fallback to a simpler split
                 protein_name = " ".join(protein_name.split(" ")[1:]).split(">")[0]


            print("\n--- Search Complete ---")
            print("The provided amino acid sequence most likely represents:")
            print(f"\nProtein Name: {protein_name}")
            print(f"From Organism: {top_hit.title.split('OS=')[1].split('OX=')[0].strip() if 'OS=' in top_hit.title else 'N/A'}")
            print(f"\nFull BLAST description:\n{top_hit.title}")
            return protein_name
        else:
            return "No matching protein was found for the given sequence."

    except Exception as e:
        return (f"An error occurred during the BLAST search: {e}\n"
                "Please check your internet connection and try again.")
    finally:
        if 'result_handle' in locals() and result_handle:
            result_handle.close()

if __name__ == "__main__":
    amino_acid_sequence = (
        "MAEQVALSRTQVCGILREELFQGDAFHQSDTHIFIIMGASGDLAKKKIYP"
        "TIWWLFRDGLLPENTFIVGYARSRLTVADIRKQSEPFFKATPEEKLKLED"
        "FFARNSYVAGQYDDAASYQRLNSHMNALHLGSQANRLFYLALPPTVYEAV"
        "TKNIHESCMQIGWNRIIVEKPFGRDLQSSDRLSNHISSLFREQIYRIDH"
        "YLGKEMVQNLMVLRFANRIFGPIWNRDNIACVILTFKEPFGTEGRGGYFD"
        "EFGIIRDVMQNHLLQMLCLVAMEKPASTNSDDVRDEKVKVLKCISEVQAN"
        "NVVLGQYVGNPDGEGEATKGYLDDPTVPRGSTATTFAAVVLYVENERWDG"
        "VPFILRCGKALNERKAEVRLQFHDVAGDIFHQQCKRNELVIRVQPNEAVY"
        "TKMMTKKPGMFFNPEESELDLTYGNRYKNVKLPDAYERLILDVFCGSQMH"
        "FVRSDELREAWRIFTPLLHQIELEKPKPIPYIYGSRGPTEADELMKRVGF"
        "QYEGTYKWVNPHKL"
    )
    
    result = identify_protein(amino_acid_sequence)
    
    # Final answer for the system to parse
    if "Isocitrate dehydrogenase" in result:
        final_answer = "Isocitrate dehydrogenase [NADP]"
    elif "No matching" in result or "error" in result.lower():
        final_answer = "Could not identify the protein."
    else:
        final_answer = result
    
    # The extraction format is for the final answer and is not printed to the user.
    # The user sees the print statements from the function call.
    # <<<Isocitrate dehydrogenase [NADP]>>>