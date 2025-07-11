import requests
import time
import re
import json

def identify_protein_sequence(raw_sequence):
    """
    Identifies a protein sequence using the NCBI BLAST API.

    Args:
        raw_sequence (str): A string containing the amino acid sequence, possibly with spaces and newlines.
    
    Returns:
        str: The name of the top protein hit, or an error message.
    """
    # Step 1: Clean the sequence
    sequence = "".join(raw_sequence.strip().split())

    base_url = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
    
    # Step 2: Submit BLASTp request to the API
    put_payload = {
        'CMD': 'Put',
        'PROGRAM': 'blastp',
        'DATABASE': 'nr',  # Non-redundant protein database
        'QUERY': sequence,
    }
    
    print("Submitting protein sequence to NCBI BLAST...")
    
    try:
        response = requests.post(base_url, data=put_payload)
        response.raise_for_status() # Raise an exception for bad status codes (4xx or 5xx)

        # The response is a formatted HTML page, we need to extract the Request ID (RID)
        match = re.search(r'RID = (\w+)', response.text)
        if not match:
            return "Error: Could not find Request ID (RID) in BLAST submission response."
            
        rid = match.group(1)
        print(f"Search successfully submitted. Request ID: {rid}")

        # Step 3: Poll for results until the search is ready
        print("Waiting for results. This may take a minute...")
        while True:
            time.sleep(10) # Wait 10 seconds between checks
            
            status_payload = {
                'CMD': 'Get',
                'RID': rid,
                'FORMAT_OBJECT': 'SearchInfo'
            }
            status_response = requests.post(base_url, data=status_payload)
            status_response.raise_for_status()
            
            if 'Status=READY' in status_response.text:
                print("Search is ready.")
                break
            elif 'Status=WAITING' in status_response.text or 'Status=RUNNING' in status_response.text:
                print("...still searching...")
                continue
            else:
                return "Error: Unknown status or error encountered during search."

        # Step 4: Retrieve the final results in JSON format
        print("Retrieving results...")
        results_payload = {
            'CMD': 'Get',
            'RID': rid,
            'FORMAT_TYPE': 'JSON2' # Request results in JSON format
        }
        results_response = requests.post(base_url, data=results_payload)
        results_response.raise_for_status()
        
        results_json = results_response.json()

        # Step 5: Parse the JSON to find the top hit
        hits = results_json.get('BlastOutput2', [{}])[0].get('report', {}).get('results', {}).get('search', {}).get('hits', [])
        
        if not hits:
            return "No hits found for the given sequence."
        
        # The first hit in the list is the best match
        top_hit_title = hits[0]['description'][0]['title']
        return top_hit_title

    except requests.exceptions.RequestException as e:
        return f"An error occurred with the network request: {e}"
    except (KeyError, IndexError) as e:
        return f"An error occurred while parsing the BLAST results: {e}"

# The amino acid sequence provided by the user
input_sequence = """
MAEQVALSRT QVCGILREEL FQGDAFHQSD THIFIIMGAS GDLAKKKIYP
TIWWLFRDGL LPENTFIVGY ARSRLTVADI RKQSEPFFKA TPEEKLKLED
FFARNSYVAG QYDDAASYQR LNSHMNALHL GSQANRLFYL ALPPTVYEAV
TKNIHESCMQ IGWNRIIVEK PFGRDLQSSD RLSNHISSLF REDQIYRIDH
YLGKEMVQNL MVLRFANRIF GPIWNRDNIA CVILTFKEPF GTEGRGGYFD
EFGIIRDVMQ NHLLQMLCLV AMEKPASTNS DDVRDEKVKV LKCISEVQAN
NVVLGQYVGN PDGEGEATKG YLDDPTVPRG STTATFAAVV LYVENERWDG
VPFILRCGKA LNERKAEVRL QFHDVAGDIF HQQCKRNELV IRVQPNEAVY
TKMMTKKPGM FFNPEESELD LTYGNRYKNV KLPDAYERLI LDVFCGSQMH
FVRSDELREA WRIFTPLLHQ IELEKPKPIP YIYGSRGPTE ADELMKRVGF
QYEGTYKWVN PHKL
"""

# Run the identification process and print the result
protein_name = identify_protein_sequence(input_sequence)

print("\n" + "="*50)
print("           Protein Identification Result")
print("="*50)
print(f"The most likely protein for the sequence is:\n{protein_name}")
print("="*50)
<<<phosphoglycerate mutase 2 [Homo sapiens]>>>