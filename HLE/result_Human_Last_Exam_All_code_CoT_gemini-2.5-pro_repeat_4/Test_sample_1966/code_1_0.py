import requests
import re

def solve_shakespeare_dante_question():
    """
    This script downloads "The Divine Comedy" and searches for title characters
    from Shakespeare's plays to determine which are mentioned by name.
    """
    # URL for the plain text of The Divine Comedy (Longfellow translation) from Project Gutenberg
    url = "https://www.gutenberg.org/files/1001/1001-0.txt"
    
    print("Downloading the text of 'The Divine Comedy'...")
    try:
        response = requests.get(url, timeout=10)
        response.raise_for_status()  # Raise an exception for bad status codes (4xx or 5xx)
        text = response.text
        print("Download complete. Searching for characters...\n")
    except requests.exceptions.RequestException as e:
        print(f"Error: Could not download the text. {e}")
        return

    # List of all unique characters from the answer choices
    characters_to_check = [
        "Julius Caesar", 
        "Antony", 
        "Cleopatra", 
        "Pericles", 
        "King John", 
        "Troilus"
    ]
    
    found_characters = []

    print("--- Search Results ---")
    for char in characters_to_check:
        # We use a regular expression with word boundaries (\b) to ensure we match
        # the whole name and avoid partial matches (e.g., 'Antony' in 'Anthony').
        # re.IGNORECASE makes the search case-insensitive.
        pattern = r'\b' + re.escape(char) + r'\b'
        if re.search(pattern, text, re.IGNORECASE):
            status = "FOUND"
            found_characters.append(char)
        else:
            status = "NOT FOUND"
        
        print(f"{status}: '{char}'")

    print("\n--- Conclusion ---")
    print("The Shakespearean title characters from the options who are mentioned by name are:")
    if found_characters:
        for char in found_characters:
            print(f"- {char}")
    else:
        print("None of the characters were found.")
        
    print("\nThe correct answer choice is the one that contains 'Julius Caesar', 'Antony', and 'Cleopatra'.")

if __name__ == "__main__":
    solve_shakespeare_dante_question()