import requests
import xml.etree.ElementTree as ET
import time

def check_pubmed_for_peptide(peptide):
    """
    Queries the PubMed database for the number of publications
    mentioning a specific peptide in the context of integrin binding.
    """
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    # The search term looks for the peptide and the word "integrin"
    search_term = f'"{peptide}" AND "integrin"'
    params = {
        "db": "pubmed",
        "term": search_term,
        "retmode": "xml",
        "tool": "my_integrin_checker", # Be a good citizen of the web
        "email": "user@example.com"  # Let NCBI know who is asking
    }
    
    try:
        response = requests.get(base_url, params=params)
        response.raise_for_status()  # Raise an exception for bad status codes
        
        root = ET.fromstring(response.content)
        count_element = root.find("Count")
        if count_element is not None:
            return int(count_element.text)
        else:
            return 0
    except requests.exceptions.RequestException as e:
        print(f"An error occurred: {e}")
        return -1 # Indicate error
    except ET.ParseError as e:
        print(f"Error parsing XML: {e}")
        return -1 # Indicate error

def main():
    """
    Main function to check all peptides and print the results.
    """
    peptides = {
        'A': 'RGDMAA',
        'B': 'RGDSPSS',
        'C': 'RGDLTTP',
        'D': 'RGDQVSK',
        'E': 'RGDARGG'
    }

    print("Querying PubMed for evidence of integrin binding for each peptide...")
    print("-" * 60)
    
    results = {}
    for choice, sequence in peptides.items():
        count = check_pubmed_for_peptide(sequence)
        results[choice] = count
        print(f"Peptide {choice} ({sequence}): Found {count} relevant publication(s).")
        # NCBI recommends not making more than 3 requests per second
        time.sleep(0.5)

    print("-" * 60)
    
    # Determine the most likely candidate
    if not results or max(results.values()) == 0:
        print("\nNo definitive evidence found via this automated search.")
        print("Based on established biochemical knowledge, the peptide RGDSPSS is derived from fibronectin's cell-binding domain and is a well-characterized integrin binder.")
        best_choice = 'B'
    else:
        # Find the choice with the maximum count
        best_choice = max(results, key=results.get)
        print(f"\nConclusion: Peptide {best_choice} ({peptides[best_choice]}) has the most supporting evidence in published literature.")

if __name__ == "__main__":
    main()