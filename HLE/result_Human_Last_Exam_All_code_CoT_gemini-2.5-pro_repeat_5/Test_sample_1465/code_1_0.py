import requests
import re
import time

def predict_oligomeric_state(sequence: str) -> str:
    """
    Submits a sequence to the LOGICOIL web server and parses the result
    to find the predicted oligomeric state.
    """
    url = "http://coiledcoils.chm.bris.ac.uk/logicoil/cgi-bin/logicoil.cgi"
    payload = {
        'seq': sequence,
        'parallel': 'par',  # Assume parallel coiled-coils, which is most common
        'outform': 'long'
    }
    
    try:
        # It's good practice to set a timeout for web requests
        response = requests.post(url, data=payload, timeout=30)
        # Raise an exception if the request was unsuccessful
        response.raise_for_status()
        
        # Use a regular expression to find the prediction in the HTML response
        # The key phrase is "The best supported assignment is that of a parallel..."
        match = re.search(r"The best supported assignment is that of a parallel (\w+)", response.text)
        
        if match:
            state = match.group(1).lower()
            if state == "dimer":
                return "2"
            elif state == "trimer":
                return "3"
            elif state == "tetramer":
                return "4"
            else:
                return f"Unknown ({state})"
        else:
            return "Prediction not found"
            
    except requests.exceptions.RequestException as e:
        # Handle potential network errors or timeouts
        print(f"An error occurred while processing sequence '{sequence[:10]}...': {e}")
        return "Error"

def main():
    """
    Main function to process the list of sequences and print the results.
    """
    sequences = [
        "EIAQALKEIAKALKEIAWALKEIAQALK",
        "EIAALKQEIAALKKENAALKQEIAALKQ",
        "EIAAIKQEIAAIKKEIAAIKWEIAAIKQ",
        "EIQKQLKEIQKQLKEIQWQLKEIQKQLK",
        "EIAQTLKEIAKTLKEIAWTLKEIAQTLK",
    ]
    
    print("Querying LOGICOIL server for oligomeric state predictions...")
    
    results = []
    for i, seq in enumerate(sequences):
        # We add a small delay to be respectful to the web server
        if i > 0:
            time.sleep(1)
        
        result = predict_oligomeric_state(seq)
        results.append(result)
        
    print("\nPrediction complete.")
    # The final output includes each number as requested
    print(f"The predicted oligomeric states are: {', '.join(results)}")

if __name__ == "__main__":
    main()