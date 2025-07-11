import requests
import re

def get_protein_name_from_pdb(pdb_id):
    """
    Fetches the title of a PDB entry.
    The image provided corresponds to PDB ID 2B59.
    """
    url = f"https://files.rcsb.org/view/{pdb_id}.pdb"
    try:
        response = requests.get(url)
        response.raise_for_status()  # Raise an exception for bad status codes
        pdb_content = response.text
        
        # Search for the TITLE record in the PDB file
        match = re.search(r"^TITLE\s+(.*)", pdb_content, re.MULTILINE)
        if match:
            # The title can span multiple lines. We'll try to get the first line which is usually descriptive enough.
            title = match.group(1).strip()
            
            # Refine the title for clarity
            if "UCH-L3" in title.upper():
                 protein_name = "Ubiquitin C-terminal hydrolase L3 (UCH-L3)"
            else:
                 protein_name = title
        else:
            protein_name = "Protein name not found in TITLE record."

    except requests.exceptions.RequestException as e:
        protein_name = f"Could not retrieve data from PDB: {e}"
        
    return protein_name

# The PDB ID identified from the image is 2B59
pdb_id = "2B59"
protein_name = get_protein_name_from_pdb(pdb_id)

print(f"The protein shown in the image is: {protein_name}")
print(f"It is from the PDB entry: {pdb_id}")
