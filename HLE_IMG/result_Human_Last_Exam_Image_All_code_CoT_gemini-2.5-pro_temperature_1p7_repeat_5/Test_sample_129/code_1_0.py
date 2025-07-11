# First, you might need to install biopython:
# pip install biopython

import os
from Bio.PDB import PDBList, PDBParser

def get_protein_names_from_pdb(pdb_id):
    """
    Downloads a PDB file, extracts the names of the molecules within it,
    and then cleans up the downloaded file.
    
    Args:
        pdb_id (str): The 4-character PDB identifier.
    """
    pdbl = PDBList()
    pdb_filename = ""
    try:
        # Download the PDB file, store it in the current directory
        print(f"Fetching data for PDB ID: {pdb_id}...")
        pdb_filename = pdbl.retrieve_pdb_file(pdb_id, pdir='.', file_format='pdb')
        
        # Check if the file was downloaded successfully
        if not os.path.exists(pdb_filename):
            print(f"Error: PDB file '{pdb_id}' could not be downloaded.")
            return

        # Parse the PDB file
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure(pdb_id, pdb_filename)
        
        # Extract compound information from the header
        compound_info = structure.header.get('compound', {})
        
        print("\n--- Protein Information ---")
        if compound_info:
            for mol_id, details in compound_info.items():
                molecule_name = details.get('molecule', 'N/A')
                chain = details.get('chain', 'N/A')
                print(f"Molecule Name: {molecule_name}")
                print(f"Associated Chain(s): {chain}")
        else:
            print("No compound information could be extracted from the PDB file header.")

    except Exception as e:
        print(f"An error occurred: {e}")
    finally:
        # Clean up by removing the downloaded PDB file
        if pdb_filename and os.path.exists(pdb_filename):
            os.remove(pdb_filename)
            # Also remove the obsolete entry if it exists
            obsolete_file = os.path.join('.', f"obsolete.dat")
            if os.path.exists(obsolete_file):
                 os.remove(obsolete_file)


if __name__ == '__main__':
    # The PDB ID corresponding to the image is 1A2Z.
    target_pdb_id = "1A2Z"
    get_protein_names_from_pdb(target_pdb_id)