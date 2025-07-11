import re
import requests

def find_unrhymed_chaucer_word():
    """
    This script finds which word from a given list does not appear in
    Chaucer's "The Book of the Duchess", and therefore cannot be rhymed.
    """
    # 1. Fetch the text of "The Book of the Duchess".
    # We will use the Project Gutenberg plain text file.
    url = "https://www.gutenberg.org/cache/epub/2245/pg2245.txt"
    try:
        response = requests.get(url)
        response.raise_for_status()  # This will raise an error for a bad response
        full_text = response.text
    except requests.exceptions.RequestException as e:
        print(f"Could not retrieve the poem from the internet. Error: {e}")
        return

    # 2. Isolate just the text of "The Book of the Duchess" from the file.
    # The file contains multiple works. We'll find the start and end markers.
    try:
        start_marker = "THE BOOK OF THE DUCHESS"
        end_marker = "Explicit the Boke of the Duchesse."
        start_index = full_text.index(start_marker)
        end_index = full_text.index(end_marker, start_index)
        poem_text = full_text[start_index:end_index]
    except ValueError:
        print("Could not isolate 'The Book of the Duchess' from the text file.")
        print("Analysis will run on the full downloaded text, which might be less accurate.")
        poem_text = full_text
    
    # Use a lowercase version for case-insensitive search
    poem_text_lower = poem_text.lower()

    # 3. Define the words to check from the answer choices.
    choices = {
        'A': 'wente',
        'B': 'here',
        'C': 'fool',
        'D': 'hool',
        'E': 'countour'
    }

    print("Checking for the presence of each word in 'The Book of the Duchess'...\n")
    
    not_found_choice = None

    # 4. Search for each word in the text.
    for key, word in choices.items():
        # We use \b for word boundaries to match the whole word only.
        pattern = r'\b' + re.escape(word) + r'\b'
        if re.search(pattern, poem_text_lower):
            print(f"Result for '{word}' (Choice {key}): Found in the text.")
        else:
            print(f"Result for '{word}' (Choice {key}): NOT found in the text.")
            not_found_choice = key
    
    print("\n-------------------------------------------------------------")
    if not_found_choice:
        print(f"Conclusion: The word '{choices[not_found_choice]}' (Choice {not_found_choice}) does not appear in the poem.")
        print("Therefore, Chaucer could not have made a rhyme with it.")
    else:
        print("Conclusion: All words were found. A more complex analysis would be needed.")

find_unrhymed_chaucer_word()
<<<E>>>