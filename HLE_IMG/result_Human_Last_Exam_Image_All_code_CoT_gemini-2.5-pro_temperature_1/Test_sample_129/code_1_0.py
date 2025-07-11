import os
from Bio.PDB import PDBList, PDBParser

def get_protein_name_from_pdb(pdb_id):
    """
    Downloads and parses a PDB file to find the protein's name.

    Args:
        pdb_id (str): The 4-character PDB ID.

    Returns:
        str: The name of the protein, or an error message.
    """
    pdbl = PDBList()
    pdb_filename = ""
    try:
        # Download the PDB file in PDB format
        pdb_filename = pdbl.retrieve_pdb_file(pdb_id, pdir='.', file_format='pdb')

        # Check if the download was successful
        if not os.path.exists(pdb_filename):
            return f"Failed to download PDB file for {pdb_id}"

        # Parse the PDB file
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure(pdb_id, pdb_filename)

        # Extract the name from the header
        protein_name = "Unknown"
        if 'name' in structure.header:
            protein_name = structure.header['name'].strip().capitalize()
        elif 'head' in structure.header:
            protein_name = structure.header['head'].strip().capitalize()
            
        return protein_name

    finally:
        # Clean up by removing the downloaded file
        if pdb_filename and os.path.exists(pdb_filename):
            os.remove(pdb_filename)
        # BioPython might create an 'obsolete' directory
        obsolete_dir = os.path.join('.', 'obsolete')
        if os.path.isdir(obsolete_dir):
            try:
                os.rmdir(obsolete_dir)
            except OSError:
                pass # Directory might not be empty, we can ignore

# The PDB ID for the protein in the image is 1UBQ.
protein_pdb_id = "1UBQ"
protein_name = get_protein_name_from_pdb(protein_pdb_id)

print(f"The protein with PDB ID '{protein_pdb_id}' is named: {protein_name}")
