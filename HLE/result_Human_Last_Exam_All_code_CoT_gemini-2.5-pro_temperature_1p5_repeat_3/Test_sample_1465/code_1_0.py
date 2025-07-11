import requests
import re
import time

def predict_oligomeric_state(sequence, seq_num):
    """
    Submits a sequence to the LOGICOIL web server and parses the result.
    
    Args:
        sequence (str): The protein sequence to analyze.
        seq_num (int): The number of the sequence for printing purposes.

    Returns:
        int or str: The predicted oligomeric state as a number (2, 3, 4) or an error message.
    """
    # URL and data payload for the LOGICOIL web server
    url = "http://coiledcoils.chm.bris.ac.uk/logicoil/cgi-bin/logicoil.cgi"
    payload = {
        'seq': sequence,
        'win': '28',
        'out': 'normal',
        'prot': 'parallel',
        'mat': 'new'
    }
    
    state_map = {"DIMER": 2, "TRIMER": 3, "TETRAMER": 4}

    try:
        # We add a small delay to be respectful to the server
        time.sleep(1)
        
        print(f"Analyzing Sequence {seq_num}: {sequence}")
        response = requests.post(url, data=payload, timeout=60)
        response.raise_for_status()  # This will raise an error for bad responses (4xx or 5xx)

        # The server returns the prediction in the format: "...predicted to be a <strong>STATE</strong>"
        # We use a regular expression to find and extract the state.
        match = re.search(r"predicted to be a <strong>(\w+)<\/strong>", response.text)

        if match:
            state_text = match.group(1).upper()
            if state_text in state_map:
                state_number = state_map[state_text]
                print(f"Result for Sequence {seq_num}: The predicted state is {state_text}, which corresponds to the number {state_number}.\n")
                return state_number
            else:
                return f"Unknown state '{state_text}' found."
        else:
            return "Could not find prediction in the server response."
            
    except requests.exceptions.RequestException as e:
        return f"An error occurred: {e}"

def main():
    """
    Main function to process all sequences and print the final result.
    """
    # List of protein sequences from the user
    sequences = [
        "EIAQALKEIAKALKEIAWALKEIAQALK",
        "EIAALKQEIAALKKENAALKQEIAALKQ",
        "EIAAIKQEIAAIKKEIAAIKWEIAAIKQ",
        "EIQKQLKEIQKQLKEIQWQLKEIQKQLK",
        "EIAQTLKEIAKTLKEIAWTLKEIAQTLK",
    ]
    
    # Predict the state for each sequence
    predicted_states = []
    for i, seq in enumerate(sequences):
        state = predict_oligomeric_state(seq, i + 1)
        if isinstance(state, int):
            predicted_states.append(state)
        else:
            print(f"Could not determine state for sequence {i+1}: {state}")
            predicted_states.append("Error")

    print("--------------------------------------------------")
    print("The final predicted oligomeric states are:")
    
    # Final output as requested
    final_output = ','.join(map(str, predicted_states))
    print(final_output)

if __name__ == '__main__':
    main()
