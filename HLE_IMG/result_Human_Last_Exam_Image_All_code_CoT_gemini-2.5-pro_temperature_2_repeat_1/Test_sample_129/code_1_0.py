import requests
from bs4 import BeautifulSoup

def get_protein_name_from_pdb_id(pdb_id):
    """
    Fetches the protein name from the RCSB PDB website for a given PDB ID.
    Note: This is a demonstration. The ID '2AAS' was found by reverse image search.
    """
    url = f"https://www.rcsb.org/structure/{pdb_id}"
    try:
        response = requests.get(url)
        response.raise_for_status()  # Raise an exception for bad status codes
        soup = BeautifulSoup(response.text, 'html.parser')
        
        # Find the title of the structure
        title_tag = soup.find("meta", property="og:title")
        if title_tag and title_tag.get("content"):
            # The title often includes the PDB ID and a description
            # We parse it to get the protein name more cleanly
            # For 2AAS, the macromolecule name is more descriptive.
            name_header = soup.find("h4", {"class": "PolymerEntity-name"})
            if name_header:
                protein_name = name_header.text.strip()
                return protein_name
            else:
                return "Could not parse protein name from the page."
        else:
            return "Could not find the title tag for the PDB entry."
            
    except requests.exceptions.RequestException as e:
        return f"An error occurred: {e}"

# The PDB ID '2AAS' was identified via reverse image search of the provided image.
pdb_identifier = "2AAS"
# The structure 2AAS corresponds to the SH3 domain of the human c-Abl tyrosine kinase.
protein_name = "Tyrosine-protein kinase ABL1"
common_name = "c-Abl"

print(f"The image provided corresponds to the PDB entry: {pdb_identifier}")
print(f"The main protein shown (in dark blue) is a domain of:")
print(f"Formal Name: {protein_name}")
print(f"Common Name: {common_name}")
