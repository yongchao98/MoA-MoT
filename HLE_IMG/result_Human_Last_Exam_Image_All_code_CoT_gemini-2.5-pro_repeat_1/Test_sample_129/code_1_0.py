try:
    from Bio.PDB import PDBList, PDBParser
except ImportError:
    print("Biopython is not installed. Please install it using: pip install biopython")
    exit()

import os

def get_protein_info(pdb_id):
    """
    Fetches a PDB file and prints its official name and compound details.
    """
    print(f"Attempting to fetch information for PDB ID: {pdb_id}...")
    
    pdbl = PDBList()
    # retrieve_pdb_file will download the file and return its path
    # It will be named pdb{pdb_id}.ent and placed in the specified directory
    pdb_file_path = pdbl.retrieve_pdb_file(pdb_id, pdir='.', file_format='pdb')

    if not os.path.exists(pdb_file_path):
        print(f"Error: Could not download PDB file for ID '{pdb_id}'.")
        return

    parser = PDBParser(QUIET=True)
    try:
        structure = parser.get_structure(pdb_id, pdb_file_path)
        
        # The 'name' field in the header provides the title of the PDB entry
        protein_name = structure.header.get('name', 'Name not found in header.')
        print("\n--- Protein Information ---")
        print(f"Official Name: {protein_name}")
        
        # The 'compound' field gives details about each molecule in the structure
        compound_info = structure.header.get('compound', {})
        if compound_info:
            print("\nComponent Molecules:")
            for key, value in compound_info.items():
                mol_name = value.get('molecule', 'Unknown molecule')
                chain_id = value.get('chain', 'N/A')
                print(f"- Molecule {key}: {mol_name} (Chain: {chain_id})")
        print("-------------------------\n")

    except Exception as e:
        print(f"An error occurred while parsing the PDB file: {e}")
    finally:
        # Clean up by removing the downloaded PDB file
        if os.path.exists(pdb_file_path):
            os.remove(pdb_file_path)
            print(f"Cleaned up downloaded file: {os.path.basename(pdb_file_path)}")

# The PDB ID for the protein complex in the image is 1FXT.
protein_pdb_id = "1FXT"
get_protein_info(protein_pdb_id)