import requests
from bs4 import BeautifulSoup

def get_protein_name_from_pdb(pdb_id):
    """
    Fetches the protein name from the RCSB PDB website for a given PDB ID.
    """
    url = f"https://www.rcsb.org/structure/{pdb_id}"
    print(f"Attempting to fetch data from: {url}")
    
    try:
        # Use a user-agent to mimic a browser
        headers = {
            'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36'
        }
        response = requests.get(url, headers=headers, timeout=15)
        response.raise_for_status()  # Raise an exception for bad status codes

        soup = BeautifulSoup(response.text, 'html.parser')

        # Find the h4 tag with the specific ID for the structure title
        title_element = soup.find('h4', id='structure_title')

        if title_element:
            protein_name = title_element.get_text(strip=True)
            return protein_name
        else:
            # Fallback: try to find by class if ID fails
            title_element = soup.find('h4', class_='text-uppercase')
            if title_element:
                 protein_name = title_element.get_text(strip=True)
                 return protein_name
            return f"Could not find the title element for PDB ID {pdb_id}."

    except requests.exceptions.RequestException as e:
        return f"An error occurred during the web request: {e}"
    except Exception as e:
        return f"An unexpected error occurred: {e}"

# The PDB ID for the protein in the image is 2B5A.
pdb_id = "2B5A"
protein_name = get_protein_name_from_pdb(pdb_id)

print("-" * 30)
print(f"The PDB ID for the structure is: {pdb_id}")
print(f"The name of the structure is: {protein_name}")
print("-" * 30)