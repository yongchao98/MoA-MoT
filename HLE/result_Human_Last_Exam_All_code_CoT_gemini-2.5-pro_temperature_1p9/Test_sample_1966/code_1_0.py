import requests
import re

def find_shakespearean_characters_in_dante():
    """
    Downloads The Divine Comedy and searches for names of Shakespearean
    title characters mentioned in the multiple-choice options.
    """
    # URL for the plain text of The Divine Comedy (Translated by H. F. Cary)
    url = "https://www.gutenberg.org/files/8800/8800-0.txt"
    
    print(f"Downloading The Divine Comedy from {url}...")
    
    try:
        response = requests.get(url, timeout=10)
        # Raise an exception for HTTP errors (e.g., 404, 500)
        response.raise_for_status()
        text = response.text
        print("Download complete. Searching for characters...\n")
    except requests.exceptions.RequestException as e:
        print(f"Error: Could not download the text. {e}")
        return

    # Characters from the answer choices. We'll search for their most common form.
    # Note: Dante refers to Julius Caesar simply as "Caesar".
    characters_to_search = {
        "Julius Caesar": "Caesar",
        "Cleopatra": "Cleopatra",
        "Troilus": "Troilus",
        "Antony": "Antony",
        "Pericles": "Pericles",
        "King John": "John" # Search for "John" as King John is unlikely
    }

    found_characters = []

    for full_name, search_term in characters_to_search.items():
        # Use regex to find whole words, ignoring case, to avoid partial matches like 'Antony' in 'Anthony'.
        # \b ensures we are matching a whole word.
        if re.search(r'\b' + re.escape(search_term) + r'\b', text, re.IGNORECASE):
            print(f"- {full_name}: Found")
            found_characters.append(full_name)
        else:
            print(f"- {full_name}: Not Found")
            
    print("\nBased on the text analysis, the title characters from Shakespeare's plays")
    print("who are mentioned by name in The Divine Comedy are:")
    if found_characters:
        # Join the found characters for a clean output sentence
        result_string = " and ".join(found_characters)
        print(result_string)
    else:
        print("None of the characters from the list were found.")

if __name__ == "__main__":
    find_shakespearean_characters_in_dante()
<<<D>>>