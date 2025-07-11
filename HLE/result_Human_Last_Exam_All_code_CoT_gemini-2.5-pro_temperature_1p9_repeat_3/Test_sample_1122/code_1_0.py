import urllib.request
import urllib.parse
import time
import re
import xml.etree.ElementTree as ET

def find_protein_identity(sequence):
    """
    Finds the identity of a protein sequence using NCBI BLAST.
    """
    print("Starting protein identification process...")
    
    # Step 1: Submit the BLASTp job to NCBI
    print("Step 1: Submitting sequence to NCBI BLAST...")
    blast_url = 'https://blast.ncbi.nlm.nih.gov/Blast.cgi'
    params = {
        'CMD': 'Put',
        'PROGRAM': 'blastp',
        'DATABASE': 'nr',
        'QUERY': sequence
    }
    
    try:
        data = urllib.parse.urlencode(params).encode('ascii')
        req = urllib.request.Request(blast_url, data)
        with urllib.request.urlopen(req) as response:
            response_html = response.read().decode('utf-8')
    except Exception as e:
        print(f"Error submitting job to NCBI: {e}")
        return

    # Step 2: Extract the Request ID (RID) from the response
    rid_match = re.search(r'RID = (\w+)', response_html)
    if not rid_match:
        print("Could not retrieve a Request ID from NCBI. Aborting.")
        return
    rid = rid_match.group(1)
    print(f"   -> Success! Job submitted with Request ID: {rid}")

    # Step 3: Poll NCBI to check if the job is finished
    print("Step 2: Waiting for the search to complete (this may take a minute)...")
    while True:
        time.sleep(10)  # Wait before polling to be considerate to NCBI servers
        
        try:
            status_url = f"{blast_url}?CMD=Get&RID={rid}&FORMAT_OBJECT=SearchInfo"
            with urllib.request.urlopen(status_url) as response:
                status_info = response.read().decode('utf-8')
        except Exception as e:
            print(f"Error checking job status: {e}")
            time.sleep(20) # Wait longer if there's an error
            continue

        status_match = re.search(r'Status=(\w+)', status_info)
        if status_match:
            status = status_match.group(1)
            if status == "READY":
                print("   -> Search is complete.")
                break
            elif status == "FAILED":
                print("   -> Search failed at NCBI. Aborting.")
                return
            else: # Status is WAITING or UNKNOWN
                print("   -> Still waiting...")
        else:
            print("   -> Could not determine job status, retrying...")

    # Step 4: Fetch the final results in XML format
    print("Step 3: Fetching and parsing the results...")
    try:
        results_url = f"{blast_url}?CMD=Get&RID={rid}&FORMAT_TYPE=XML"
        with urllib.request.urlopen(results_url) as response:
            xml_results = response.read().decode('utf-8')
    except Exception as e:
        print(f"Error fetching final results: {e}")
        return

    # Step 5: Parse the XML and print the top hit
    try:
        root = ET.fromstring(xml_results)
        # The path to the top hit's description in the BLAST XML
        hit_def_element = root.find('.//Hit_def')
        
        if hit_def_element is not None and hit_def_element.text:
            protein_name = hit_def_element.text
            print("\n--- Identification Complete ---")
            print("The most likely identity of the protein is:")
            print(protein_name)
        else:
            print("Could not find a protein match in the results.")
            
    except ET.ParseError:
        print("Failed to parse XML results from NCBI.")

if __name__ == '__main__':
    # The amino acid sequence provided by the user
    amino_acid_sequence = (
        "MAEQVALSRTQVCGILREELFQGDAFHQSDTHIFIIMGASGDLAKKKIYP"
        "TIWWLFRDGLLPENTFIVGYARSRLTVADIRKQSEPFFKATPEEKLKLED"
        "FFARNSYVAGQYDDAASYQRLNSHMNALHLGSQANRLFYLALPPTVYEAV"
        "TKNIHESCMQIGWNRIIVEKPFGRDLQSSDRLSNHISSLFREDQIYRIDH"
        "YLGKEMVQNLMVLRFANRIFGPIWNRDNIACVILTFKEPFGTEGRGGYFD"
        "EFGIIRDVMQNHLLQMLCLVAMEKPASTNSDDVRDEKVKVLKCISEVQAN"
        "NVVLGQYVGNPDGEGEATKGYLDDPTVPRGSTTATFAAVVLYVENERWDG"
        "VPFILRCGKALNERKAEVRLQFHDVAGDIFHQQCKRNELVIRVQPNEAVY"
        "TKMMTKKPGMFFNPEESELDLTYGNRYKNVKLPDAYERLILDVFCGSQMH"
        "FVRSDELREAWRIFTPLLHQIELEKPKPIPYIYGSRGPTEADELMKRVGF"
        "QYEGTYKWVNPHKL"
    )
    find_protein_identity(amino_acid_sequence)