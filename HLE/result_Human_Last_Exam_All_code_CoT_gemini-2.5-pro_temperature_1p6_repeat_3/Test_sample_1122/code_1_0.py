import requests
import time
import re
import xml.etree.ElementTree as ET

def identify_protein_sequence():
    """
    Identifies a protein by querying its amino acid sequence against the NCBI BLAST database.
    """
    # Step 1: Define and clean the amino acid sequence
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
    sequence = "".join(raw_sequence.split())

    print(f"Querying for protein with sequence starting with: {sequence[:30]}...")
    print("-" * 70)

    # Step 2: Use the NCBI BLAST API
    base_url = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"

    # --- Stage 1: SUBMIT search and get Request ID (RID) ---
    print("1. Submitting BLASTP search to NCBI...")
    submit_params = {
        'CMD': 'Put',
        'PROGRAM': 'blastp',
        'DATABASE': 'nr',  # Non-Redundant Protein Sequence Database
        'SEQUENCE': sequence
    }
    try:
        response = requests.post(base_url, params=submit_params, timeout=30)
        response.raise_for_status()  # Raise an exception for bad status codes (4xx or 5xx)
    except requests.exceptions.RequestException as e:
        print(f"Error submitting BLAST request: {e}")
        return

    # Parse the response to get the RID and RTOE (estimated time)
    rid_match = re.search(r'RID = (\S+)', response.text)
    rtoe_match = re.search(r'RTOE = (\d+)', response.text)

    if not rid_match:
        print("Could not find Request ID (RID) in the NCBI response.")
        print("Response Text:", response.text)
        return

    rid = rid_match.group(1)
    rtoe = int(rtoe_match.group(1)) if rtoe_match else 20
    print(f"   - Success! Request ID: {rid}")

    # --- Stage 2: POLL for results until the search is ready ---
    print(f"2. Waiting for the search to complete (estimated: {rtoe}s)...")
    time.sleep(rtoe)

    check_params = {'CMD': 'Get', 'RID': rid, 'FORMAT_OBJECT': 'SearchInfo'}
    while True:
        try:
            check_response = requests.get(base_url, params=check_params, timeout=30)
            check_response.raise_for_status()
            if "Status=READY" in check_response.text:
                print("   - Search is complete.")
                break
            elif "Status=FAILED" in check_response.text:
                print("   - Search failed. Please check the sequence or try again later.")
                return
            else:
                print("   - Search is still running, checking again in 10 seconds...")
                time.sleep(10)
        except requests.exceptions.RequestException as e:
            print(f"Error checking BLAST status: {e}")
            return

    # --- Stage 3: RETRIEVE and parse the final results ---
    print("3. Retrieving and parsing results...")
    retrieve_params = {'CMD': 'Get', 'RID': rid, 'FORMAT_TYPE': 'XML'}

    try:
        result_response = requests.get(base_url, params=retrieve_params, timeout=30)
        result_response.raise_for_status()
    except requests.exceptions.RequestException as e:
        print(f"Error retrieving BLAST results: {e}")
        return

    # Parse the XML to find the name of the top hit
    try:
        root = ET.fromstring(result_response.content)
        hit = root.find('.//Hit')  # The first Hit is the best match
        if hit is not None:
            hit_def = hit.find('Hit_def').text
            hsp = hit.find('.//Hsp')
            e_value = float(hsp.find('Hsp_evalue').text)
            identity = int(hsp.find('Hsp_identity').text)
            positives = int(hsp.find('Hsp_positive').text)
            query_len = int(hsp.find('Hsp_align-len').text)
            identity_percent = (identity / query_len) * 100
            positives_percent = (positives / query_len) * 100

            print("\n" + "=" * 70)
            print("                TOP PROTEIN MATCH FOUND")
            print("=" * 70)
            print(f"Protein Name: {hit_def}")
            print(f"E-value:      {e_value:.2e} (A measure of significance; closer to 0 is better)")
            print(f"Identity:     {identity}/{query_len} ({identity_percent:.2f}%)")
            print(f"Positives:    {positives}/{query_len} ({positives_percent:.2f}%)")
            print("-" * 70)

        else:
            print("No significant matches found for the provided sequence.")

    except ET.ParseError:
        print("Could not parse the XML results from NCBI.")

if __name__ == '__main__':
    identify_protein_sequence()