# First, ensure you have biopython installed:
# pip install biopython

import os
from Bio.PDB import PDBList, PDBParser

def get_protein_name(pdb_id: str):
    """
    Downloads a PDB file and prints the name of the protein.

    Args:
        pdb_id (str): The 4-character PDB ID.
    """
    # Create a PDBList object to manage downloading files
    pdbl = PDBList()

    # Download the PDB file (in .ent format) to the current directory
    # The 'overwrite=True' flag will replace the file if it already exists.
    # The 'file_format' is set to 'pdb' which downloads a .ent file.
    file_path = pdbl.retrieve_pdb_file(pdb_id, file_format='pdb', pdir='.')
    
    # Check if the file was downloaded successfully
    if not os.path.exists(file_path):
        print(f"Failed to download PDB file for ID: {pdb_id}")
        return

    try:
        # Create a PDBParser object
        parser = PDBParser(QUIET=True)

        # Parse the structure from the downloaded file
        structure = parser.get_structure(pdb_id, file_path)

        # The protein name is stored in the header dictionary of the structure object
        protein_name = structure.header.get('name', 'Name not found')
        
        print(f"PDB ID: {pdb_id}")
        print(f"Protein Name: {protein_name}")

    except Exception as e:
        print(f"An error occurred: {e}")
    finally:
        # Clean up by removing the downloaded file
        if os.path.exists(file_path):
            os.remove(file_path)

# The PDB ID for the protein in the image is 2B5A
protein_id = "2B5A"
get_protein_name(protein_id)