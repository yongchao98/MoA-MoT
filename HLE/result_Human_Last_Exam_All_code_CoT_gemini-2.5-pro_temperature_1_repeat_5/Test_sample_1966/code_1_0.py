import requests
import re

def find_shakespearean_characters_in_dante():
    """
    Fetches 'The Divine Comedy' and searches for the names of specific
    Shakespearean title characters mentioned in the multiple-choice options.
    """
    # URL for the plain text of The Divine Comedy (H. W. Longfellow translation)
    url = "https://www.gutenberg.org/files/1001/1001-0.txt"
    
    try:
        response = requests.get(url)
        # Raise an exception if the download fails
        response.raise_for_status()
        text = response.text
    except requests.exceptions.RequestException as e:
        print(f"Could not retrieve the text from the URL. Error: {e}")
        return

    # Characters from the answer choices to search for.
    # The regex pattern \b ensures we match whole words only.
    characters = {
        'Julius Caesar': r'\bCaesar\b',
        'Pericles': r'\bPericles\b',
        'Cleopatra': r'\bCleopatra\b',
        'King John': r'\bKing John\b',
        'Troilus': r'\bTroilus\b',
        'Antony': r'\bAntony\b'
    }

    print("Searching for Shakespearean title characters in 'The Divine Comedy'...")
    print("-" * 60)

    found_characters = []

    for name, pattern in characters.items():
        # re.IGNORECASE makes the search case-insensitive
        if re.search(pattern, text, re.IGNORECASE):
            found_characters.append(name)
            print(f"[FOUND]     '{name}' is mentioned by name.")
        else:
            print(f"[NOT FOUND] '{name}' is not mentioned by name.")
            
    print("-" * 60)
    
    # Sort the list for consistent output
    found_characters.sort()
    
    # Final summary statement as required by the prompt
    print("Based on the text, the Shakespearean title characters mentioned by name are: " + ", ".join(found_characters) + ".")

find_shakespearean_characters_in_dante()
