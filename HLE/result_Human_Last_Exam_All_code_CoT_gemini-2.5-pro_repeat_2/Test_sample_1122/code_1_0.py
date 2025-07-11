import requests
import time
import sys

def identify_protein_sequence():
    """
    Identifies a protein by performing a BLAST search on its amino acid sequence
    using the EBI REST API.
    """
    # The amino acid sequence to identify
    sequence = (
        "MAEQVALSRTQVCGILREELFQGDAFHQSDTHIFIIMGASGDLAKKKIYPTIWWLFRDGLLPENTFIVGY"
        "ARSRLTVADIRKQSEPFFKATPEEKLKLEDFFARNSYVAGQYDDAASYQRLNSHMNALHLGSQANRLFYL"
        "ALPPTVYEAVTKNIHESCMQIGWNRIIVEKPFGRDLQSSDRLSNHISSLFREDOIYRIDHYLGKEMVQNL"
        "MVLRFANRIFGPIWNRDNIACVILTFKEPFGTEGRGGYFDEFGIIRDVMQNHLLQMLCLVAMEKPASTNS"
        "DDVRDEKVKVLKCISEVQANNVVLGQYVGNPDGEGEATKGYLDDPTVPRGSTATTFAAVVLYVENERWDG"
        "VPFILRCGKALNERKAEVRLQFHDVAGDIFHQQCKRNELVIRVQPNEAVYTKKMMTKKPGMFFNPEESEL"
        "DLTYGNRYKNVKLPDAYERLILDVFCGSQMHHFVRSDELREAWRIFTPLLHQIELEKPKPIPYIYGSRGP"
        "TEADELMKRVGFQYEGTYKWVNPHKL"
    )

    # EBI BLAST API endpoint details
    base_url = "https://www.ebi.ac.uk/Tools/services/rest/blastp"
    run_url = f"{base_url}/run"
    
    # It's good practice to provide a user-agent and a placeholder email
    headers = {
        'User-Agent': 'Python-Client/1.0',
    }
    payload = {
        'email': 'user@example.com',
        'database': 'uniprotkb',
        'program': 'blastp',
        'sequence': sequence
    }

    # 1. Submit the BLAST job
    try:
        print("Submitting BLAST job...")
        run_response = requests.post(run_url, headers=headers, data=payload)
        run_response.raise_for_status()
        job_id = run_response.text
        print(f"Job submitted successfully. Job ID: {job_id}")
    except requests.exceptions.RequestException as e:
        print(f"Error submitting job: {e}", file=sys.stderr)
        return

    # 2. Poll for job status
    status_url = f"{base_url}/status/{job_id}"
    print("Waiting for results... This may take a minute.")
    while True:
        try:
            status_response = requests.get(status_url, headers=headers)
            status_response.raise_for_status()
            status = status_response.text
            
            if status == 'FINISHED':
                print("Job finished.")
                break
            elif status in ('NOT_FOUND', 'ERROR', 'FAILURE'):
                print(f"Job failed with status: {status}", file=sys.stderr)
                return
            
            # Wait for 10 seconds before checking again
            time.sleep(10)
        except requests.exceptions.RequestException as e:
            print(f"Error checking status: {e}", file=sys.stderr)
            return

    # 3. Fetch the results
    result_url = f"{base_url}/result/{job_id}/out?format=json"
    print("Fetching results...")
    try:
        result_response = requests.get(result_url, headers=headers)
        result_response.raise_for_status()
        results = result_response.json()
    except (requests.exceptions.RequestException, ValueError) as e:
        print(f"Error fetching or parsing results: {e}", file=sys.stderr)
        return

    # 4. Parse and print the top hit
    try:
        if not results['results']['hits']:
            print("No significant matches found for the sequence.")
            return

        top_hit = results['results']['hits'][0]
        accession = top_hit['hit_acc']
        description = top_hit['hit_description']
        # The first HSP (High-scoring Segment Pair) usually contains the key stats
        hsp = top_hit['hit_hsps'][0]
        score = hsp['hsp_score']
        e_value = hsp['hsp_expect']
        identity = hsp['hsp_identity']
        length = hsp['hsp_align_len']
        
        print("\n--- Top BLAST Hit ---")
        print(f"Protein Name: {description}")
        print(f"Accession ID: {accession}")
        print(f"Score: {score}")
        print(f"E-value: {e_value}")
        print(f"Identity: {identity}/{length} ({100 * int(identity) / int(length):.2f}%)")

    except (IndexError, KeyError) as e:
        print(f"Could not parse results. Unexpected format. Error: {e}", file=sys.stderr)

if __name__ == '__main__':
    identify_protein_sequence()