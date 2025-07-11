import requests
import re
import time

def identify_protein_sequence():
    """
    Identifies a protein by submitting its amino acid sequence to NCBI BLAST
    and parsing the top hit from the results.
    """
    # The amino acid sequence provided by the user
    sequence = """MAEQVALSRT QVCGILREEL FQGDAFHQSD THIFIIMGAS GDLAKKKIYP 
    TIWWLFRDGL LPENTFIVGY ARSRLTVADI RKQSEPFFKA TPEEKLKLED 
    FFARNSYVAG QYDDAASYQR LNSHMNALHL GSQANRLFYL ALPPTVYEAV 
    TKNIHESCMQ IGWNRIIVEK PFGRDLQSSD RLSNHISSLF REDQIYRIDH 
    YLGKEMVQNL MVLRFANRIF GPIWNRDNIA CVILTFKEPF GTEGRGGYFD 
    EFGIIRDVMQ NHLLQMLCLV AMEKPASTNS DDVRDEKVKV LKCISEVQAN 
    NVVLGQYVGN PDGEGEATKG YLDDPTVPRG STTATFAAVV LYVENERWDG 
    VPFILRCGKA LNERKAEVRL QFHDVAGDIF HQQCKRNELV IRVQPNEAVY 
    TKMMTKKPGM FFNPEESELD LTYGNRYKNV KLPDAYERLI LDVFCGSQMH 
    FVRSDELREA WRIFTPLLHQ IELEKPKPIP YIYGSRGPTE ADELMKRVGF 
    QYEGTYKWVN PHKL"""

    # Remove newlines and whitespace to create a single continuous sequence string
    cleaned_sequence = "".join(sequence.split())

    # --- Step 1: Submit the BLAST search to NCBI API ---
    print("Submitting sequence to NCBI BLAST...")
    put_url = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
    put_payload = {
        'CMD': 'Put',
        'PROGRAM': 'blastp',
        'DATABASE': 'nr',
        'QUERY': cleaned_sequence
    }
    
    try:
        response = requests.post(put_url, data=put_payload, timeout=30)
        response.raise_for_status()
    except requests.RequestException as e:
        print(f"Error during BLAST submission: {e}")
        return

    # Extract the Request ID (RID) from the submission response
    rid_match = re.search(r"RID = (\S+)", response.text)
    if not rid_match:
        print("Error: Could not find Request ID (RID) in the BLAST response.")
        return
        
    rid = rid_match.group(1)
    print(f"Search successfully submitted. Request ID: {rid}")
    print("Waiting for search to complete...")

    # --- Step 2: Poll server until results are ready ---
    get_url = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
    while True:
        # Wait a few seconds before checking the status to avoid overloading the server
        time.sleep(10)
        
        get_payload = {'CMD': 'Get', 'RID': rid, 'FORMAT_OBJECT': 'SearchInfo'}
        try:
            status_response = requests.get(get_url, params=get_payload, timeout=30)
            status_response.raise_for_status()
        except requests.RequestException as e:
            print(f"Error checking BLAST status: {e}")
            return

        # Check the status within the response
        if "Status=READY" in status_response.text:
            if "ThereAreHits=yes" in status_response.text:
                print("Search complete. Retrieving results...")
                break
            else:
                print("Search complete, but no matching proteins were found.")
                return
        elif "Status=FAILED" in status_response.text:
            print("Search failed. This may be due to an invalid sequence or a server error.")
            return
        # If status is still "WAITING", the loop will continue

    # --- Step 3: Retrieve and parse the final results ---
    results_payload = {'CMD': 'Get', 'RID': rid, 'FORMAT_TYPE': 'Text'}
    try:
        results_response = requests.get(get_url, params=results_payload, timeout=60)
        results_response.raise_for_status()
        results_text = results_response.text
    except requests.RequestException as e:
        print(f"Error retrieving BLAST results: {e}")
        return

    # Find and parse the top hit from the summary table in the results
    try:
        summary_section = results_text.split("Sequences producing significant alignments:")[1]
        summary_section = summary_section.split("ALIGNMENTS")[0]

        # Find the first line in the summary table, which corresponds to the best match
        top_hit_line = ""
        for line in summary_section.strip().split('\n'):
            # Skip empty lines and the header line of the table
            if line.strip() and not line.strip().startswith(('(bits)', 'Description')):
                top_hit_line = line
                break

        # Extract details using a regular expression
        match = re.match(r'(.+?)\s+([\d\.]+)\s+([\de\.\-]+)', top_hit_line)
        if match:
            description = match.group(1).strip()
            score = match.group(2)
            e_value = match.group(3)
            
            print("\n--- Top Hit Identified ---")
            print(f"Description: {description}")
            print(f"Score = {score}")
            print(f"E-value = {e_value}")
        else:
            print("\nCould not parse top hit details, showing raw line:")
            print(top_hit_line)

    except (IndexError, AttributeError):
        print("Error: Failed to parse the BLAST results format.")

if __name__ == "__main__":
    identify_protein_sequence()