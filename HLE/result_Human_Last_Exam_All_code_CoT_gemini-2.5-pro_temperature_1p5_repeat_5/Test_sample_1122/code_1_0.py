import requests
import re
import time
import json

def identify_protein_sequence(sequence):
    """
    Identifies a protein sequence using the NCBI BLAST API.
    
    Args:
        sequence (str): The amino acid sequence to identify.
    """
    
    # Step 1: Submit the sequence to the BLAST API
    print("Submitting sequence to NCBI BLAST...")
    put_url = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
    put_params = {
        'CMD': 'Put',
        'PROGRAM': 'blastp',
        'DATABASE': 'nr', # Non-redundant protein sequences
        'QUERY': sequence
    }
    
    try:
        response = requests.post(put_url, data=put_params)
        response.raise_for_status()
    except requests.exceptions.RequestException as e:
        print(f"Error submitting sequence to BLAST: {e}")
        return

    # Extract the Request ID (RID) from the response HTML
    rid_match = re.search(r"RID = ([\w-]+)", response.text)
    if not rid_match:
        print("Could not find BLAST Request ID (RID).")
        print("Response Text:", response.text[:500]) # Print first 500 chars for debugging
        return
        
    rid = rid_match.group(1)
    print(f"Search submitted. Request ID (RID): {rid}")

    # Step 2: Poll for the results to be ready
    print("Waiting for results...")
    get_url = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
    
    while True:
        try:
            get_params = {'CMD': 'Get', 'RID': rid, 'FORMAT_OBJECT': 'SearchInfo'}
            status_response = requests.get(get_url, params=get_params)
            status_response.raise_for_status()

            # Check if the search is ready
            if "Status=READY" in status_response.text:
                print("Search complete.")
                break
            # Check for failures
            elif "Status=FAILED" in status_response.text:
                print("Search failed. Please try again later.")
                return
            
            # Wait before checking again
            time.sleep(10)

        except requests.exceptions.RequestException as e:
            print(f"Error checking BLAST status: {e}")
            return
            
    # Step 3: Fetch the final results in JSON format
    print("Fetching results...")
    try:
        result_params = {'CMD': 'Get', 'RID': rid, 'FORMAT_TYPE': 'JSON2'}
        result_response = requests.get(get_url, params=result_params)
        result_response.raise_for_status()
        
        results = result_response.json()
        
        # Parse the JSON to find the top hit
        search_results = results.get('BlastOutput2', [])
        if not search_results:
            print("No results found in the BLAST output.")
            return

        hits = search_results[0].get('report', {}).get('results', {}).get('search', {}).get('hits', [])
        if not hits:
            print("No matching proteins found.")
            return

        # The first hit in the list is the best match
        top_hit = hits[0]
        description = top_hit['description'][0]['title']
        
        print("\n--- Top Match Found ---")
        print(f"The most likely protein is: {description}")

    except (requests.exceptions.RequestException, json.JSONDecodeError) as e:
        print(f"Error fetching or parsing BLAST results: {e}")
    except (KeyError, IndexError) as e:
        print(f"Could not parse the protein name from the BLAST results. Error: {e}")
        print("Received data structure might have changed.")

if __name__ == '__main__':
    amino_acid_sequence = (
        "MAEQVALSRT QVCGILREEL FQGDAFHQSD THIFIIMGAS GDLAKKKIYP "
        "TIWWLFRDGL LPENTFIVGY ARSRLTVADI RKQSEPFFKA TPEEKLKLED "
        "FFARNSYVAG QYDDAASYQR LNSHMNALHL GSQANRLFYL ALPPTVYEAV "
        "TKNIHESCMQ IGWNRIIVEK PFGRDLQSSD RLSNHISSLF REDQIYRIDH "
        "YLGKEMVQNL MVLRFANRIF GPIWNRDNIA CVILTFKEPF GTEGRGGYFD "
        "EFGIIRDVMQ NHLLQMLCLV AMEKPASTNS DDVRDEKVKV LKCISEVQAN "
        "NVVLGQYVGN PDGEGEATKG YLDDPTVPRG STTATFAAVV LYVENERWDG "
        "VPFILRCGKA LNERKAEVRL QFHDVAGDIF HQQCKRNELV IRVQPNEAVY "
        "TKMMTKKPGM FFNPEESELD LTYGNRYKNV KLPDAYERLI LDVFCGSQMH "
        "FVRSDELREA WRIFTPLLHQ IELEKPKPIP YIYGSRGPTE ADELMKRVGF "
        "QYEGTYKWVN PHKL"
    ).replace(" ", "") # Remove spaces
    
    identify_protein_sequence(amino_acid_sequence)