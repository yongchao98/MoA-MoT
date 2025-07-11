import requests
import time

def identify_protein_sequence():
    """
    Identifies a protein sequence by performing a BLAST search against the UniProt database.
    """
    # The amino acid sequence provided by the user.
    sequence = """
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

    # 1. Clean up the sequence string by removing newlines and spaces.
    cleaned_sequence = "".join(sequence.strip().split())

    # 2. Define the UniProt REST API endpoint for submitting a BLAST job.
    run_url = "https://rest.uniprot.org/blast/run"

    # 3. Set up the parameters and submit the job via a POST request.
    params = {
        "program": "blastp",
        "database": "uniprotkb_swissprot",
        "sequence": cleaned_sequence
    }
    print("Submitting sequence to UniProt BLAST...")
    try:
        submit_req = requests.post(run_url, json=params)
        submit_req.raise_for_status() # Raise an exception for bad status codes
        job_id = submit_req.json().get("jobId")
        if not job_id:
            print("Error: Could not retrieve Job ID from UniProt.")
            print("Response:", submit_req.json())
            return
    except requests.exceptions.RequestException as e:
        print(f"Error submitting job to UniProt: {e}")
        return

    print(f"Job submitted successfully. Job ID: {job_id}")

    # 4. Poll the status endpoint until the job is finished.
    status_url = f"https://rest.uniprot.org/blast/status/{job_id}"
    while True:
        try:
            status_req = requests.get(status_url)
            status_req.raise_for_status()
            status_data = status_req.json()
            job_status = status_data.get("jobStatus")

            if job_status == "FINISHED":
                print("BLAST search completed.")
                break
            elif job_status in ["RUNNING", "QUEUED"]:
                print(f"Job is {job_status.lower()}... checking again in 5 seconds.")
                time.sleep(5)
            else:
                print(f"Job failed with status: {job_status}")
                print("Details:", status_data)
                return
        except requests.exceptions.RequestException as e:
            print(f"Error checking job status: {e}")
            time.sleep(10) # Wait longer if there's a connection issue

    # 5. Fetch and parse the final results.
    results_url = f"https://rest.uniprot.org/blast/results/{job_id}"
    print(f"Fetching results from {results_url}")
    try:
        results_req = requests.get(results_url)
        results_req.raise_for_status()
        results_data = results_req.json()
    except requests.exceptions.RequestException as e:
        print(f"Error fetching results: {e}")
        return

    # 6. Display the details of the top hit.
    if not results_data.get("results"):
        print("No significant hits were found for the provided sequence.")
    else:
        top_hit = results_data["results"][0]

        hit_accession = top_hit.get("uniProtkbAccession", "N/A")
        
        protein_info = top_hit.get("protein", {})
        full_name_info = protein_info.get("recommendedName", {}).get("fullName", {})
        protein_name = full_name_info.get("value", "N/A")

        organism_info = top_hit.get("organism", {})
        organism_name = organism_info.get("scientificName", "N/A")
        
        identity = top_hit.get("match", {}).get("sequenceIdentity", 0.0) * 100

        print("\n--- Top BLAST Hit Found ---")
        print(f"Protein Name: {protein_name}")
        print(f"Accession ID: {hit_accession}")
        print(f"Organism:     {organism_name}")
        print(f"Identity:     {identity:.2f}%")

if __name__ == '__main__':
    identify_protein_sequence()