import requests
import re
import time

def find_protein_name(sequence):
    """
    Identifies a protein by its amino acid sequence using NCBI BLASTp.
    """
    print("Submitting the amino acid sequence to NCBI BLAST for identification...")
    
    # Step 1: Submit the BLASTp job via NCBI API
    put_url = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
    payload = {
        'CMD': 'Put',
        'PROGRAM': 'blastp',
        'DATABASE': 'nr',
        'QUERY': sequence
    }
    
    try:
        response = requests.post(put_url, data=payload)
        response.raise_for_status() # Raise an exception for bad status codes
    except requests.exceptions.RequestException as e:
        print(f"An error occurred while submitting the BLAST request: {e}")
        return

    # Step 2: Extract the Request ID (RID) from the response
    rid_match = re.search(r'RID = (\S+)', response.text)
    if not rid_match:
        print("Failed to get a Request ID from NCBI. The service might be busy or the request was invalid.")
        return
        
    rid = rid_match.group(1)
    print(f"Search successfully submitted. Request ID (RID) is: {rid}")

    # Step 3: Poll for results until the search is ready
    print("Waiting for the search to complete...")
    get_url = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
    while True:
        time.sleep(5) # Wait before polling
        status_payload = {
            'CMD': 'Get',
            'FORMAT_OBJECT': 'SearchInfo',
            'RID': rid
        }
        try:
            status_response = requests.get(get_url, params=status_payload)
            status_response.raise_for_status()
            status_match = re.search(r'Status=(\w+)', status_response.text)
            
            if status_match:
                status = status_match.group(1)
                if status == "READY":
                    print("Search has completed.")
                    break
                elif status == "FAILED":
                    print("The BLAST search failed. Please try again later.")
                    return
                # If status is "WAITING", the loop continues
        except requests.exceptions.RequestException as e:
            print(f"An error occurred while checking search status: {e}")
            return

    # Step 4: Retrieve and parse the final results
    print("Retrieving and parsing results...")
    results_payload = {
        'CMD': 'Get',
        'RID': rid,
        'FORMAT_TYPE': 'Text' # Request results as plain text
    }
    try:
        results_response = requests.get(get_url, params=results_payload)
        results_response.raise_for_status()
        results_text = results_response.text
        
        # Find the description of the top hit
        # The top hit is the first line after the "Sequences producing significant alignments:" header
        if "Sequences producing significant alignments:" in results_text:
            alignments_section = results_text.split("Sequences producing significant alignments:")[1]
            lines = alignments_section.strip().split('\n')
            if lines:
                # The first line contains the name of the top hit
                top_hit_line = lines[0]
                # The name is the part before the score, which is separated by multiple spaces
                protein_name = top_hit_line.split('   ')[0].strip()
                print("\n--- Found Protein ---")
                print(f"The sequence most likely represents the protein: {protein_name}")
            else:
                print("Could not find any significant alignments in the results.")
        else:
            print("No significant similarity found.")
            
    except (requests.exceptions.RequestException, IndexError) as e:
        print(f"An error occurred while retrieving or parsing the results: {e}")


if __name__ == '__main__':
    # The amino acid sequence provided by the user
    amino_acid_sequence = (
        "MAEQVALSRTQVCGILREELFQGDAFHQSDTHIFIIMGASGDLAKKKIYPTIWWLFRDGLLPENTFIVGYARSRLTVADIRKQSEPFFKATPEEKLKLEDFFARNSYVAGQYDDAASYQRLNSHMNALHLGSQANRLFYLALPPTVYEAVTKNIHESCMQIGWNRIIVEKPFGRDLQSSDRLSNHISSLFREDQIYRIDHYLGKEMVQNLMVLRFANRIFGPIWNRDNIACVILTFKEPFGTEGRGGYFDEFGIIRDVMQNHLLQMLCLVAMEKPASTNSDDVRDEKVKVLKCISEVQANNVVLGQYVGNPIDGEGEATKGYLDDPTVPRGSTATFAAVVLYVENERWDGVPFILRCGKALNERKAEVRLQFHDVAGDIFHQQCKRNELVIRVQPNEAVYTKKMMTKKPGMFFNPEESELDLTYGNRYKNVKLPDAYERLILDVFCGSQMHFVRSDELREAWRIFTPLLHQIELEKPKPIPIYGSRGPTEADELMKRVGFQYEGTYKWVNPHKL"
    )
    find_protein_name(amino_acid_sequence)