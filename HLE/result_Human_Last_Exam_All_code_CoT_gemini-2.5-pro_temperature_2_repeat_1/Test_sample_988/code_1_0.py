import xml.etree.ElementTree as ET
from Bio import Entrez

def find_initial_antioxidant_response():
    """
    Searches the PubMed database to find which antioxidants are initially
    activated in Microcystis aeruginosa CAAT 2005-3 in response to
    high temperature exposure.
    """
    # Provide an email to NCBI. This is a requirement for using the Entrez API.
    Entrez.email = "ai.assistant@example.com"

    # Define the search terms based on the user's question
    search_query = '"Microcystis aeruginosa CAAT 2005-3" AND "high temperature" AND "oxidative stress"'

    print(f"Searching PubMed database with the query: '{search_query}'...")
    print("-" * 70)

    try:
        # Search PubMed for publication IDs (PMIDs)
        search_handle = Entrez.esearch(db="pubmed", term=search_query, retmax=1)
        search_results = Entrez.read(search_handle)
        search_handle.close()

        id_list = search_results["IdList"]

        if not id_list:
            print("No relevant publications found for this specific query.")
            return

        pmid = id_list[0]
        print(f"Found relevant publication with PMID: {pmid}. Fetching abstract...")

        # Fetch the full record for the top publication ID
        fetch_handle = Entrez.efetch(db="pubmed", id=pmid, rettype="xml", retmode="text")
        xml_data = fetch_handle.read()
        fetch_handle.close()

        # Parse the XML to extract the title and abstract
        root = ET.fromstring(xml_data)
        article = root.find(".//PubmedArticle")

        if article is not None:
            title = article.findtext(".//ArticleTitle")
            abstract = ""
            # Find the AbstractText element, which might not always exist
            abstract_element = article.find(".//AbstractText")
            if abstract_element is not None:
                abstract = abstract_element.text

            print(f"\nTitle: {title}\n")
            print(f"Abstract Snippet:\n{abstract}\n")
            print("-" * 70)
            print("Conclusion from the Abstract:")
            if "enzymatic antioxidants" in abstract.lower() and "initially triggered" in abstract.lower():
                print("The abstract states that 'enzymatic antioxidants' ... is the main antioxidant mechanism to be 'initially triggered'.")
                print("This directly corresponds to Answer Choice C.")
            else:
                print("Analysis of the abstract is required to determine the answer.")

    except Exception as e:
        print(f"An error occurred while connecting to NCBI: {e}")

# Run the search
find_initial_antioxidant_response()