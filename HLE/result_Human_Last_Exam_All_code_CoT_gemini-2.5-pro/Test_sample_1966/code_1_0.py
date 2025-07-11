import requests
import re

def find_shakespearean_characters_in_dante():
    """
    This script finds which of Shakespeare's title characters from a given list
    are mentioned in Dante's 'The Divine Comedy' and determines the correct
    answer from a set of choices.
    """
    # URL for the plain text of "The Divine Comedy" (Longfellow translation)
    url = "https://www.gutenberg.org/cache/epub/8800/pg8800.txt"
    
    # All unique characters from the answer choices
    characters_to_search = ["Julius Caesar", "Pericles", "Cleopatra", "King John", "Troilus", "Antony"]

    print("Step 1: Downloading the text of 'The Divine Comedy'...")
    try:
        response = requests.get(url)
        response.raise_for_status()  # This will raise an error for a bad response (4xx or 5xx)
        text = response.text
        print("Download complete.")
    except requests.exceptions.RequestException as e:
        print(f"Error: Could not download the text. {e}")
        return

    print("\nStep 2: Searching for each character in the text...")
    
    found_characters = []
    for character in characters_to_search:
        # For multi-word names like "Julius Caesar", a simple string search is effective.
        # For single names, use regex with word boundaries (\b) to avoid partial matches.
        # The search is case-insensitive (re.IGNORECASE).
        if len(character.split()) > 1:
            if character.lower() in text.lower():
                found_characters.append(character)
                print(f"- Found: {character}")
            else:
                print(f"- Not Found: {character}")
        else:
            if re.search(r'\b' + re.escape(character) + r'\b', text, re.IGNORECASE):
                found_characters.append(character)
                print(f"- Found: {character}")
            else:
                print(f"- Not Found: {character}")

    print("\nStep 3: Comparing the list of found characters with the answer choices...")

    # The choices provided in the problem
    choices = {
        "A": {"Julius Caesar", "Pericles"},
        "B": {"Julius Caesar", "Cleopatra", "King John"},
        "C": {"Julius Caesar", "Troilus", "Antony"},
        "D": {"Julius Caesar", "Cleopatra"},
        "E": {"Julius Caesar", "Antony", "Cleopatra"}
    }
    
    found_set = set(found_characters)
    correct_choice = None

    for key, value in choices.items():
        if value == found_set:
            correct_choice = key
            break
            
    if correct_choice:
        print(f"\nThe characters found ({', '.join(sorted(list(found_set)))}) exactly match the set in choice {correct_choice}.")
    else:
        print("\nNone of the choices exactly match the list of characters found.")


    print("\nFinal Answer:")
    if correct_choice:
        print(f"<<<{correct_choice}>>>")
    else:
        # Fallback in case of an unexpected result
        print("<<<Could not determine a matching answer.>>>")

# Run the function
find_shakespearean_characters_in_dante()