import requests
import time
import sys

def identify_protein_sequence():
    """
    Identifies a protein sequence by performing a BLAST search against the UniProt database.
    """
    # The amino acid sequence provided by the user.
    raw_sequence = """
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
    # Clean up the sequence: remove spaces and newlines.
    sequence = "".join(raw_sequence.strip().split())

    base_url = "https://www.uniprot.org/uniprot/api"
    submit_url = f"{base_url}/blast/run"
    
    # The payload for the POST request includes the program (blastp for protein),
    # the database to search against, and the sequence itself.
    payload = {
        "program": "blastp",
        "stype": "protein",
        "database": "uniprotkb_refprotswissprot", # Search against curated Swiss-Prot
        "sequence": sequence,
    }

    try:
        print("Submitting protein sequence to UniProt BLAST...")
        submit_response = requests.post(submit_url, data=payload)
        submit_response.raise_for_status()  # Raise an exception for HTTP errors
        job_id = submit_response.json()["jobId"]
        print(f"Job submitted successfully. Job ID: {job_id}")

        # Poll the status URL until the job is finished
        status_url = f"{base_url}/blast/status/{job_id}"
        while True:
            status_response = requests.get(status_url)
            status_response.raise_for_status()
            status_data = status_response.json()
            if status_data.get("jobStatus") == "FINISHED":
                print("BLAST job has finished.")
                break
            print("Job is running, checking status in 10 seconds...")
            time.sleep(10)

        # Fetch the results
        results_url = f"{base_url}/blast/results/{job_id}"
        print("Fetching BLAST results...")
        results_response = requests.get(results_url)
        results_response.raise_for_status()
        results_data = results_response.json()

        # Extract and display information for the top hit
        if results_data.get("results") and len(results_data["results"]) > 0:
            top_hit = results_data["results"][0]
            entry = top_hit.get("entry", {})
            protein_info = entry.get("proteinDescription", {})
            recommended_name = protein_info.get("recommendedName", {})
            full_name = recommended_name.get("fullName", {})
            protein_name = full_name.get("value", "N/A")

            accession = entry.get("primaryAccession", "N/A")
            identity_score = top_hit.get('sequenceIdentity', 'N/A')
            e_value = top_hit.get('evalue', 'N/A')

            print("\n--- Top Match Found ---")
            print(f"Protein Name: {protein_name}")
            print(f"UniProt Accession ID: {accession}")
            print(f"Identity: {identity_score}%")
            print(f"E-value: {e_value}")
            print("-----------------------")
        else:
            print("No significant matches found for the provided sequence.")

    except requests.exceptions.RequestException as e:
        print(f"An error occurred while communicating with the UniProt API: {e}", file=sys.stderr)
    except (KeyError, IndexError) as e:
        print(f"Could not parse the API response. It might have an unexpected format. Error: {e}", file=sys.stderr)
    except Exception as e:
        print(f"An unexpected error occurred: {e}", file=sys.stderr)

if __name__ == '__main__':
    identify_protein_sequence()